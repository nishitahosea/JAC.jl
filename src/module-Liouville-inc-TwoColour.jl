# module-Liouville-inc-TwoColour.jl

using ..Basics, ..Defaults, ..Pulse, ..PhotoExcitation, ..PhotoEmission

# Call this in getDipoleFromPhotoExcitation

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
#function getDipoleFromPhotoExcitation(level_i::Level, level_j::Level, grid::Radial.Grid)
#    # Determine final (higher energy) and initial (lower energy)
#    if level_j.energy > level_i.energy
#        final_level = level_j
#        initial_level = level_i
#    else
#        final_level = level_i
#        initial_level = level_j
#    end
#
#    omega = abs(level_j.energy - level_i.energy)
#
#    if omega < 1e-6
#        return 0.0 + 0.0im
#    end
#
#    # Convert AngularJ64 to Float64 for calculations
#    Ji = Float64(initial_level.J)
#    Jf = Float64(final_level.J)
#    delta_J = abs(Jf - Ji)
#    parity_change = (final_level.parity != initial_level.parity)
#
#    println("\n=== Dipole Debug ===")
#    println("  initial_level: J=$Ji, parity=$(initial_level.parity), energy=$(initial_level.energy)")
#    println("  final_level:   J=$Jf, parity=$(final_level.parity), energy=$(final_level.energy)")
#    println("  omega = $omega a.u.")
#    println("  delta_J = $delta_J, parity_change = $parity_change")
#    println("  E1 allowed? $(delta_J <= 1.0 && !(delta_J == 0.0 && Ji == 0.0) && parity_change)")
#
#    # Check E1 selection rules
#    if delta_J > 1.0 || (delta_J == 0.0 && Ji == 0.0) || !parity_change
#        println("  → Transition not allowed by E1 selection rules")
#        return 0.0 + 0.0im
#    end
#
#    # Use PhotoExcitation for bound-bound transitions
#    settings = PhotoExcitation.Settings(
#        [E1], [UseCoulomb], false, false, false, false,
#        LineSelection(), 0.0, 0.0, 1.0e6, Basics.ExpStokes()
#    )
#
#    channels = PhotoExcitation.determineChannels(final_level, initial_level, settings)
#    println("  Number of channels found: $(length(channels))")
#
#    if isempty(channels)
#        println("  → No channels found")
#        return 0.0 + 0.0im
#    end
#
#    for ch in channels
#        println("  Channel: multipole=$(ch.multipole), gauge=$(ch.gauge)")
#    end
#
#    line = PhotoExcitation.Line(initial_level, final_level, omega,
#                                EmProperty(0.,0.), EmProperty(0.,0.),
#                                TensorComp[], true, channels)
#    computed_line = PhotoExcitation.computeAmplitudesProperties(line, grid, settings, printout=false)
#
#    for channel in computed_line.channels
#        if channel.multipole == E1 && channel.gauge == Basics.Coulomb
#            println("  → Dipole amplitude = $(channel.amplitude)")
#            return channel.amplitude
#        end
#    end
#
#    println("  → No E1 Coulomb channel found")
#    return 0.0 + 0.0im
#end
function getDipoleFromPhotoExcitation(level_i::Level, level_j::Level, grid::Radial.Grid)
    # Determine final (higher energy) and initial (lower energy)
    if level_j.energy > level_i.energy
        final_level = level_j
        initial_level = level_i
    else
        final_level = level_i
        initial_level = level_j
    end

    omega = abs(level_j.energy - level_i.energy)
    if omega < 1e-6
        return 0.0 + 0.0im
    end

    # Use PhotoExcitation to compute the line
    settings = PhotoExcitation.Settings(
        [E1], [UseCoulomb], false, false, false, false,
        LineSelection(), 0.0, 0.0, 1.0e6, Basics.ExpStokes()
    )
    channels = PhotoExcitation.determineChannels(final_level, initial_level, settings)
    if isempty(channels)
        return 0.0 + 0.0im
    end
    line = PhotoExcitation.Line(initial_level, final_level, omega,
                                EmProperty(0.,0.), EmProperty(0.,0.),
                                TensorComp[], true, channels)
    computed_line = PhotoExcitation.computeAmplitudesProperties(line, grid, settings, printout=false)

    # Extract oscillator strength (Coulomb gauge) – gives correct magnitude
    f = computed_line.oscStrength.Coulomb
    if f <= 0.0
        return 0.0 + 0.0im
    end

    # Magnitude from oscillator strength: |d| = sqrt(3f/(2ω))
    d_mag = sqrt(3 * f / (2 * omega))

    # Extract phase from the raw amplitude of the first E1 Coulomb channel
    phase = 1.0 + 0.0im
    for channel in computed_line.channels
        if channel.multipole == E1 && channel.gauge == Basics.Coulomb
            amp = channel.amplitude
            if abs(amp) > 1e-12
                phase = amp / abs(amp)   # unit complex number (e^{iφ})
            end
            break
        end
    end

    return d_mag * phase
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
            push!(field_funcs, t -> pulse.A0 * envelope(pulse, t) * cos(pulse.omega * (t - pulse.timeDelay)))
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

    # DEBUG: Print field at t=0.1
    println("\n=== Field Debug ===")
    println("Total E(0.1) = $(total_field_func(0.1)) a.u.")

    # Build coupling matrix
    for i in 1:noLevels
        for j in 1:noLevels
            if i == j
                couplingHM[i,j] = t -> 0.0 + 0.0im
            else
                # Skip if either level is the loss channel (energy == 0.0)
                if levels[i].level.energy == 0.0 || levels[j].level.energy == 0.0
                    couplingHM[i,j] = t -> 0.0 + 0.0im
                else
                    # Ensure Hermiticity: always compute dipole from lower to higher energy
                    if levels[i].level.energy < levels[j].level.energy
                        dipole = getDipoleFromPhotoExcitation(levels[i].level, levels[j].level, grid)
                    else
                        dipole = conj(getDipoleFromPhotoExcitation(levels[j].level, levels[i].level, grid))
                    end
                    couplingHM[i,j] = t -> -dipole * total_field_func(t)
                end
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
    t = 0.0  # XUV peaks at t=0

    totalHamiltonian = [t -> couplingHM[i,j](t) + atomicHM[i,j] for i in 1:noLevels, j in 1:noLevels]
    totalH = [f(t) for f in totalHamiltonian]

    println(stream, " ")
    println(stream, "  Two-Color Hamiltonian matrix (atomic + coupling), evaluated for t=$t:")
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
            # Change from %8.6f to %10.8f to see more decimal places
            sa = sa * @sprintf("%10.8f %+10.8fim  ", real(z), imag(z))
        end
        println(sa)
    end
    println(stream, " ")
    return nothing
end


"""
`Liouville.testDipole()`
    ... test function to verify dipole calculation between 1s and 2p levels.
"""
function testDipole()
    println("\n=== Testing Dipole Calculation ===")

    # Setup
    nm = Nuclear.Model(2.0)
    grid = Radial.Grid(true)
    refConfigs = [Configuration("1s"), Configuration("2p")]
    asfSettings = AsfSettings()

    # Compute multiplet
    multiplet = SelfConsistent.performSCF(refConfigs, nm, grid, asfSettings)

    # Get 1s and 2p levels
    level1s = nothing
    level2p = nothing

    for level in multiplet.levels
        if occursin("1s", string(level.configuration))
            level1s = level
            println("Found 1s: J=$(level.J), parity=$(level.parity), energy=$(level.energy)")
        elseif occursin("2p", string(level.configuration))
            level2p = level
            println("Found 2p: J=$(level.J), parity=$(level.parity), energy=$(level.energy)")
        end
    end

    if level1s !== nothing && level2p !== nothing
        dipole = getDipoleFromPhotoExcitation(level1s, level2p, grid)
        println("Dipole moment (1s → 2p): $dipole a.u.")
    else
        println("Could not find 1s and 2p levels")
    end
end


function testField()
    println("\n=== Testing Field Magnitude ===")

    xuv_pulse = Pulse.FelPulse("GaussianSimplified", 21.0, 1e12, 10.0, 0.0, 0.1)
    nir_pulse = Pulse.FelPulse("GaussianSimplified", 1.55, 1e14, 10.0, 10.0, 0.1)

    xuv_comp = Pulse.convertPulse(xuv_pulse)
    nir_comp = Pulse.convertPulse(nir_pulse)

    println("XUV peak at t = $(xuv_comp.timeDelay) a.u.")
    println("NIR peak at t = $(nir_comp.timeDelay) a.u.")

    for t in [0.0, 100.0, 200.0, 300.0, 400.0, 413.4, 500.0]
        E_xuv = xuv_comp.A0 * envelope(xuv_comp, t) * cos(xuv_comp.omega * (t - xuv_comp.timeDelay))
        E_nir = nir_comp.A0 * envelope(nir_comp, t) * cos(nir_comp.omega * (t - nir_comp.timeDelay))
        println("t = $t: E_xuv = $E_xuv, E_nir = $E_nir, total = $(E_xuv + E_nir)")
    end
end

