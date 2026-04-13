# module-Liouville-inc-TwoColour.jl

using ..Basics, ..Defaults, ..Pulse, ..PhotoExcitation

function inspect_J(J)
    println("Type of J: ", typeof(J))
    println("Fieldnames: ", fieldnames(typeof(J)))
    for field in fieldnames(typeof(J))
        println("  $field = ", getfield(J, field))
    end
    # Try all conversion methods
    try println("  J.value = ", J.value) catch; println("  .value not available") end
    try println("  Int(J) = ", Int(J)) catch; println("  Int(J) not available") end
    try println("  Float64(J) = ", Float64(J)) catch; println("  Float64(J) not available") end
end

# Call this in getDipoleFromPhotoExcitation
inspect_J(initial_level.J)
# Define missing envelope function
function envelope(pulse::Pulse.GaussianSimplified, t::Float64)
    sigma = pulse.fwhm / (2 * sqrt(2 * log(2)))
    wa = (t - pulse.timeDelay)^2 / (2 * sigma^2)
    return exp(-wa) / (sigma * sqrt(2pi))
end

function gaussianEnvelope(t::Float64, sigma::Float64)
    wa = t^2 / (2*sigma)^2
    return exp(-wa) / (sigma * sqrt(2pi))
end

"""
`struct Liouville.TwoColourLevel`
    ... defines a struct to comprise the level information for a Liouville time evolution in the two-color scheme.
"""
struct TwoColourLevel
    leadingConfig      ::Configuration
    leadingNotation    ::String
    level              ::Level
end

function TwoColourLevel()
    TwoColourLevel(Configuration(), "xx", Level())
end

function Base.show(io::IO, level::TwoColourLevel)
    println(io, "leadingConfig:          $(level.leadingConfig)")
    println(io, "leadingNotation:        $(level.leadingNotation)")
    println(io, "level:                  $(level.level)")
end

"""
`Liouville.initializeLevels(scheme::TwoColourScheme, multiplet::Multiplet)`
    ... initialize the levels for two-color computation.
"""
function initializeLevels(scheme::TwoColourScheme, multiplet::Multiplet)
    liouvilleLevels = TwoColourLevel[];
    noLevels = length(scheme.levelSelection.indices)

    # Check proper level notations
    if length(scheme.levelNotations) - 1 != noLevels
        error("Expect $noLevels strings for level notations, got $(scheme.levelNotations)")
    end

    for (idx, index) in enumerate(scheme.levelSelection.indices)
        for level in multiplet.levels
            if index == level.index
                leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
                liouvLevel = TwoColourLevel(leadingConf, scheme.levelNotations[idx], level)
                push!(liouvilleLevels, liouvLevel)
            end
        end
    end

    # Add a loss channel (ionization continuum)
    push!(liouvilleLevels, TwoColourLevel(Configuration("[He]"), scheme.levelNotations[end], Level()))

    # Display levels
    println(" ")
    println("  Selected Two-Color levels:")
    println(" ")
    for (idx, level) in enumerate(liouvilleLevels)
        sa = "       " * string(idx) * ")  "
        sa = sa * string(level.leadingConfig) * "   "
        sa = sa * string(level.leadingNotation) * "   "
        println(sa)
    end
    println(" ")

    return liouvilleLevels
end

"""
`Liouville.initializeDensityMatrix(levels::Array{TwoColourLevel,1})`
    ... initialize the density matrix (ground state populated).
"""
function initializeDensityMatrix(levels::Array{TwoColourLevel,1})
    noLevels = length(levels)
    energies = Float64[]
    densityM = zeros(ComplexF64, noLevels, noLevels)

    for level in levels
        push!(energies, level.level.energy)
    end
    lowestEn = minimum(energies)
    idx = findfirst(==(lowestEn), energies)
    densityM[idx, idx] = 1.0

    return densityM
end

"""
`Liouville.displayDensityMatrix(stream, levels::Array{TwoColourLevel,1}, densityM::Matrix{ComplexF64})`
    ... Display the current density matrix.
"""
function displayDensityMatrix(stream, levels::Array{TwoColourLevel,1}, densityM::Matrix{ComplexF64})
    println(stream, " ")
    println(stream, "  Selected Two-Color levels and current density matrix:")
    println(stream, " ")
    for (idx, level) in enumerate(levels)
        sa = "       " * string(idx) * ")  "
        sa = sa * string(level.leadingConfig) * "   "
        sa = sa * string(level.leadingNotation) * "                          "
        # Only truncate if string is long enough
        if length(sa) >= 70
            sa = sa[1:70]
        end
        row = densityM[idx, :]
        for z in row
            sa = sa * @sprintf("%8.3f %+8.3fim  ", real(z), imag(z))
        end
        println(sa)
    end
    println(stream, " ")
    return nothing
end


"""
`Liouville.initializeAtomicHamiltonianMatrix(scheme::TwoColourScheme, levels::Array{TwoColourLevel,1})`
    ... initialize the atomic Hamiltonian matrix (diagonal energies).
"""
function initializeAtomicHamiltonianMatrix(scheme::TwoColourScheme, levels::Array{TwoColourLevel,1})
    noLevels = length(levels)
    energies = Float64[]
    hamiltonian = zeros(ComplexF64, noLevels, noLevels)

    # #=Get=# energies for all levels
    for level in levels
        push!(energies, level.level.energy)
    end

    # Find lowest energy for reference
    lowestEn = minimum(energies)

    # Set diagonal elements as excitation energies
    for n in 1:noLevels
        hamiltonian[n, n] = energies[n] - lowestEn
    end

    return hamiltonian
end



"""
`Liouville.getDipoleFromPhotoExcitation(level_i::Level, level_j::Level, grid::Radial.Grid)`
    ... computes electric dipole matrix element using JAC's PhotoExcitation module.
"""
function getDipoleFromPhotoExcitation(level_i::Level, level_j::Level, grid::Radial.Grid)
    # Determine final (higher energy) and initial (lower energy)
    if level_j.energy > level_i.energy
        final_level = level_j
        initial_level = level_i
    else
        final_level = level_i
        initial_level = level_j
    end

    println("\n=== Debug: Inspecting J values ===")
    println("initial_level.J:")
    inspect_J(initial_level.J)
    println("final_level.J:")
    inspect_J(final_level.J)

    omega = abs(level_j.energy - level_i.energy)

    if omega < 1e-6
        return 0.0 + 0.0im
    end

    # Check E1 selection rules
    # ΔJ = 0, ±1 (but not 0→0)
    # Parity must change
    delta_J = abs(final_level.J - initial_level.J)
    parity_change = (final_level.parity != initial_level.parity)

    # Convert AngularJ64 to Float64 for comparison
    delta_J_float = Float64(delta_J)
    initial_J_float = Float64(initial_level.J)

    if delta_J_float > 1.0 || (delta_J_float == 0.0 && initial_J_float == 0.0) || !parity_change
        return 0.0 + 0.0im  # Transition not allowed
    end

    try
        # PhotoEmission.amplitude expects (multipole, gauge, omega, finalLevel, initialLevel, grid)
        dipole = PhotoEmission.amplitude(E1, Basics.Coulomb, omega, final_level, initial_level, grid, printout=false)
        return dipole
    catch e
        @warn "Dipole calculation failed: $e"
        return 0.0 + 0.0im
    end
end


"""
`Liouville.initializeCouplingHamiltonianMatrix(scheme::TwoColourScheme, levels::Array{TwoColourLevel,1}, pulses::Array{Pulse.AbstractPulse, 1})`
    ... initialize the coupling Hamiltonian matrix for two-color XUV+NIR interaction.
"""
function initializeCouplingHamiltonianMatrix(scheme::TwoColourScheme, levels::Array{TwoColourLevel,1},
                                             pulses::Array{Pulse.AbstractPulse, 1}, grid::Radial.Grid)
    noLevels = length(levels)
    couplingHM = Array{Function, 2}(undef, noLevels, noLevels)

    # Create field functions for each pulse
    field_funcs = []
    for pulse in pulses
        if typeof(pulse) == Pulse.GaussianSimplified
            push!(field_funcs, t -> pulse.A0 * envelope(pulse, t) * cos(pulse.omega * (t - pulse.timeDelay) + pulse.cep))
        else
            error("Unknown pulse type: $(typeof(pulse))")
        end
    end

    total_field_func = t -> begin
        sum = 0.0
        for f in field_funcs
            sum += f(t)
        end
        return sum
    end

    # Build coupling matrix
    for i in 1:noLevels
        for j in 1:noLevels
            if i == j
                couplingHM[i,j] = t -> 0.0 + 0.0im
            else
                # Get dipole using PhotoExcitation (no grid argument)
                dipole = getDipoleFromPhotoExcitation(levels[i].level, levels[j].level, grid)
                couplingHM[i,j] = t -> -dipole * total_field_func(t)
            end
        end
    end

    return couplingHM
end


function computeInteractionMatrix(scheme::TwoColourScheme, levels::Array{TwoColourLevel,1})
    noLevels = length(levels)
    couplingHM = Array{Function, 2}(undef, noLevels, noLevels)

    # Initialize with zero functions
    for i = 1:noLevels, j = 1:noLevels
        couplingHM[i,j] = t -> 0.0 + 0.0im
    end

    for i = 1:noLevels, j = 1:noLevels
        if i == j continue end

        # 1. Calculate the static dipole matrix element <i|d|j>
        # You would use JAC's InteractionMatrix tools here
        dipoleElement = PhotoIonization.amplitude("multipole: E1", levels[i].level, levels[j].level)

        # 2. Assign the time-dependent coupling
        # Pulse 1 (XUV) and Pulse 2 (IR)
        p1 = scheme.pulses[1]
        p2 = scheme.pulses[2]

        couplingHM[i,j] = t -> dipoleElement * (Pulse.getAmplitude(t, p1) + Pulse.getAmplitude(t, p2))
    end
    return couplingHM
end


"""
`Liouville.perform(scheme::TwoColourScheme, computation::Computation; output::Bool=true)`
    ... perform two-color XUV+NIR Liouville computation.
"""
function perform(scheme::TwoColourScheme, computation::Computation; output::Bool=true)
    if output    results = Dict{String, Any}()    else    results = nothing    end

    println("")
    printstyled("Liouville.perform(): Two-Color XUV+NIR computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------ \n", color=:light_green)

    # Convert pulses
    pulses = Pulse.AbstractPulse[]
    for pulse in computation.pulses
        if typeof(pulse) == Pulse.GaussianSimplified
            push!(pulses, pulse)
        elseif typeof(pulse) == Pulse.FelPulse
            push!(pulses, Pulse.convertPulse(pulse))
        else
            error("Unknown pulse = $pulse")
        end
    end

    # Atomic structure
    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel,
                                          computation.grid, computation.asfSettings)

    # Initialize levels and density matrix
    levels = initializeLevels(scheme, multiplet)
    densityM = initializeDensityMatrix(levels)

    # ========== NEW CODE: Initialize Hamiltonians ==========
    atomicHM = initializeAtomicHamiltonianMatrix(scheme, levels)
    couplingHM = initializeCouplingHamiltonianMatrix(scheme, levels, pulses, computation.grid) # Pass the 'pulses' array here
    # ======================================================

    if computation.settings.printBefore
        displayDensityMatrix(stdout, levels, densityM)
        # ========== NEW CODE: Display Hamiltonians ==========
        displayGenericHamiltonian(stdout, levels, atomicHM, couplingHM)
        # ====================================================
    end

    println("\n  Two-Color time evolution not yet implemented.")

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
        results["pulses"] = computation.pulses
        # ========== NEW CODE: Add Hamiltonians to results ==========
        results["atomic_hamiltonian"] = atomicHM
        results["coupling_hamiltonian"] = couplingHM
        # ============================================================
    end

    println("\n> Two-Color computation setup complete ...")

    return results
end


"""
`Liouville.displayGenericHamiltonian(stream, levels::Array{TwoColourLevel,1}, atomicHM::Matrix{ComplexF64},
                                     couplingHM::Array{Function, 2})`
    ... Display the Hamiltonian matrices at a sample time t=0.1.
"""
function displayGenericHamiltonian(stream, levels::Array{TwoColourLevel,1}, atomicHM::Matrix{ComplexF64},
                                   couplingHM::Array{Function, 2})
    noLevels = length(levels)
    totalHamiltonian = [t -> couplingHM[i,j](t) + atomicHM[i,j] for i in 1:noLevels, j in 1:noLevels]

    t = 0.1
    totalH = [f(t) for f in totalHamiltonian]

    println(stream, " ")
    println(stream, "  Two-Color Hamiltonian matrix (atomic + coupling), evaluated for t=0.1:")
    println(stream, " ")
    for (idx, level) in enumerate(levels)
        sa = "       " * string(idx) * ")  "
        sa = sa * string(level.leadingConfig) * "   "
        sa = sa * string(level.leadingNotation) * "                          "
        if length(sa) >= 70
            sa = sa[1:70]
        end
        row = totalH[idx, :]
        for z in row
            sa = sa * @sprintf("%8.3f %+8.3fim  ", real(z), imag(z))
        end
        println(sa)
    end
    println(stream, " ")
    return nothing
end
