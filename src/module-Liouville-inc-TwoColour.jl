# module-Liouville-inc-TwoColour.jl

using ..Basics, ..Defaults, ..Pulse, ..PhotoExcitation, ..PhotoEmission, ..PhotoIonization, ..Continuum, ..Nuclear

# Define envelope function
function envelope(pulse::Pulse.GaussianSimplified, t::Float64)
    sigma = pulse.fwhm / (2 * sqrt(2 * log(2)))
    wa = (t - pulse.timeDelay)^2 / (2 * sigma^2)
    return exp(-wa)
end


function carrier( pulse::Pulse.GaussianSimplified, t::Float64 )
    return cos( pulse.omega * t )
end

function pulse( pulse::Pulse.GaussianSimplified, t::Float64 )
    return envelope( pulse, t ) * carrier( pulse * t )
end


"""
`struct Liouville.TwoColourLevel`
    ... defines a struct to comprise the level information for a Liouville time evolution in the two-color scheme.
"""
struct TwoColourLevel
    leadingConfig      ::Configuration
    leadingNotation    ::String
    level              ::Level
    isContinuum        ::Bool
end

function TwoColourLevel()
    TwoColourLevel(Configuration(), "xx", Level(), false)
end

function Base.show(io::IO, level::TwoColourLevel)
    println(io, "leadingConfig:          $(level.leadingConfig)")
    println(io, "leadingNotation:        $(level.leadingNotation)")
    println(io, "level:                  $(level.level)")
    println(io, "isContinuum:            $(level.isContinuum)")
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
                liouvLevel = TwoColourLevel(leadingConf, scheme.levelNotations[idx], level, false)
                push!(liouvilleLevels, liouvLevel)
            end
        end
    end

    # Add a loss channel (ionization continuum)
    push!(liouvilleLevels, TwoColourLevel(Configuration("[He]"), scheme.levelNotations[end], Level(), false))

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
    getTransitionAmplitude(level_i::Level, level_j::Level, grid::Radial.Grid, multipole::EmMultipole=E1, gauge::EmGauge=Basics.Coulomb)

Returns the complex transition amplitude (reduced matrix element) for a bound‑bound
transition between level_i and level_j for the specified multipole and gauge.
"""
function getTransitionAmplitude(level_i::Level, level_j::Level, grid::Radial.Grid,
                                multipole=E1, gauge=Basics.Coulomb)
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

    # Settings with the requested multipole and gauge
    settings = PhotoExcitation.Settings(
        [multipole], [gauge], false, false, false, false,
        LineSelection(), 0.0, 0.0, 1.0e6, Basics.ExpStokes()
    )

    # Determine allowed channels
    channels = PhotoExcitation.determineChannels(final_level, initial_level, settings)
    if isempty(channels)
        return 0.0 + 0.0im
    end

    # Create and compute the line
    line = PhotoExcitation.Line(initial_level, final_level, omega,
                                EmProperty(0.,0.), EmProperty(0.,0.),
                                TensorComp[], true, channels)
    computed_line = PhotoExcitation.computeAmplitudesProperties(line, grid, settings, printout=false)

    # Extract the amplitude - use the exact gauge from the computed channel
    for channel in computed_line.channels
        if channel.multipole == multipole && string(channel.gauge) == "Coulomb"
            # The amplitude from PhotoExcitation is the reduced matrix element.
            # For E1 transitions, this is the dipole matrix element in a.u.
            return channel.amplitude
        end
    end

    return 0.0 + 0.0im
end

function getDipoleFromPhotoExcitation(level_i::Level, level_j::Level, grid::Radial.Grid)
    return getTransitionAmplitude(level_i, level_j, grid, E1, Basics.Coulomb)
end


function getDipoleBoundContinuum(boundLevel::Level, contLevel::Level, grid::Radial.Grid, nm::Nuclear.Model)
    # Energy difference (photon energy)
    omega = abs(contLevel.energy - boundLevel.energy)
    if omega < 1e-8
        return 0.0 + 0.0im
    end

    # Extract kappa from the continuum level (assuming it contains exactly one extra electron)
    # This depends on how you build the continuum level; adapt as needed.
    extraElectrons = Basics.extraElectrons(contLevel.basis, boundLevel.basis)  # you may need a helper
    if length(extraElectrons) != 1
        error("Continuum level does not have exactly one extra electron.")
    end
    kappa = extraElectrons[1].subshell.kappa

    # Build the photoionization channel
    sym = LevelSymmetry(contLevel.J, contLevel.parity)
    channel = PhotoIonization.Channel(E1, Basics.Coulomb, kappa, sym, 0.0, 0.0+0.0im)

    # Compute amplitude
    amp = PhotoIonization.amplitude("photoionization", channel, omega, contLevel, boundLevel, grid)
    return amp
end

"""
`Liouville.initializeCouplingHamiltonianMatrix(scheme::TwoColourScheme, levels::Array{TwoColourLevel,1}, pulses::Array{Pulse.AbstractPulse, 1})`
    ... initialize the coupling Hamiltonian matrix for two-color XUV+NIR interaction.
"""
function initializeCouplingHamiltonianMatrix(scheme::TwoColourScheme, levels::Array{TwoColourLevel,1},
                                             pulses::Array{Pulse.AbstractPulse,1}, grid::Radial.Grid, nm::Nuclear.Model)
    noLevels = length(levels)
    couplingHM = Array{Function,2}(undef, noLevels, noLevels)

    # Create field functions for each pulse (sum of all fields)
    field_funcs = []
    for pulse in pulses
        if typeof(pulse) == Pulse.GaussianSimplified
            # E(t) = A0 * envelope(t) * cos(ω(t - t_d))
            push!(field_funcs, t -> pulse.A0 * envelope(pulse, t) * cos(pulse.omega * (t - pulse.timeDelay)))
        else
            error("Unknown pulse type: $(typeof(pulse))")
        end
    end

    total_field_func = t -> sum(f(t) for f in field_funcs)

    # Build coupling matrix – only bound‑bound couplings
    for i in 1:noLevels, j in 1:noLevels
        if i == j
            couplingHM[i,j] = t -> 0.0 + 0.0im
            continue
        end

        li = levels[i]
        lj = levels[j]

        # Skip the loss channel (energy == 0.0)
        if li.level.energy == 0.0 || lj.level.energy == 0.0
            couplingHM[i,j] = t -> 0.0 + 0.0im
            continue
        end

        # Only bound‑bound couplings
        # Ensure Hermiticity: always compute dipole from lower to higher energy
        if li.level.energy < lj.level.energy
            dipole = getDipoleFromPhotoExcitation(li.level, lj.level, grid)
        else
            dipole = conj(getDipoleFromPhotoExcitation(lj.level, li.level, grid))
        end

        # V_ij = -d_ij * E_total(t)
        couplingHM[i,j] = t -> -dipole * total_field_func(t)
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
`Liouville.perform(scheme::Liouville.StimulatedTwoColour, computation::Liouville.Computation; output::Bool=true)`
    ... to perform a Liouville time-evolution computation for a Two Colour scheme. For output=true, a dictionary
        is returned from which the relevant results can be can easily accessed by proper keys.
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
    noLevels = length(levels)
    densityM = initializeDensityMatrix(levels)

    # Initialize Hamiltonians
    atomicHM = initializeAtomicHamiltonianMatrix(scheme, levels)
    couplingHM = initializeCouplingHamiltonianMatrix(scheme, levels, pulses, computation.grid, computation.nuclearModel)

    if computation.settings.printBefore
        displayDensityMatrix(stdout, levels, densityM)
        displayGenericHamiltonian(stdout, levels, atomicHM, couplingHM)
    end

    # Extract NIR pulse
    nir_pulse = nothing
    for pulse in pulses
        if pulse.omega < 0.1
            nir_pulse = pulse
            break
        end
    end

    if nir_pulse === nothing
        error("No NIR pulse found in pulses")
    end

    nir_omega = nir_pulse.omega

    # Pre‑compute peak ionization rates
    gamma_peak = zeros(Float64, noLevels)
    for i in 1:noLevels
        if levels[i].level.energy < 0.0
            gamma_peak[i] = get_total_ionization_rate(levels[i].level, nir_omega, computation.nuclearModel, computation.grid)
        end
    end

    # Store results (no ODE solving here)
    if output
        results["levels"] = levels
        results["initial_density"] = densityM
        results["pulses"] = computation.pulses
        results["atomic_hamiltonian"] = atomicHM
        results["coupling_hamiltonian"] = couplingHM
        results["gamma_peak"] = gamma_peak
        results["nir_pulse"] = nir_pulse
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

function evaluateHamiltonianMatrix(t, atomicHM, couplingHM)
    n = size(atomicHM, 1)
    H = copy(atomicHM)
    for i in 1:n, j in 1:n
        H[i,j] += couplingHM[i,j](t)
    end
    return H
end


function get_total_ionization_rate(level::Level, omega::Float64, nm::Nuclear.Model, grid::Radial.Grid)
    Ji = level.J
    parity_i = level.parity

    # Convert Ji to Float64 for arithmetic
    Ji_float = Float64(Ji)

    # Choose final J such that ΔJ = 0, ±1 (but not 0→0)
    if Ji_float == 0.0
        Jf = AngularJ64(1, 1)   # J = 1
    elseif Ji_float == 0.5
        Jf = AngularJ64(1, 2)   # J = 1/2 (ΔJ=0)
    else
        Jf = AngularJ64(Int(2*Ji_float + 2), 2)   # J = Ji + 1
    end

    # Flip parity correctly
    if parity_i == Basics.plus
        parity_f = Basics.minus
    else
        parity_f = Basics.plus
    end

    dummy_level = Level(Jf, AngularM64(0), parity_f, 0, 0.0, 0.0, false, Basis(), Float64[])
    dummy_multiplet = Multiplet("dummy", [dummy_level])
    initial_multiplet = Multiplet("initial", [level])

    settings = PhotoIonization.Settings(
        [E1], [UseCoulomb], [omega], Float64[], Float64[], Float64[],
        false, false, false, false, false, false, LineSelection(),
        Basics.ExpStokes(), 0.0, [0,1,2,3,4,5]
    )

    lines = PhotoIonization.determineLines(dummy_multiplet, initial_multiplet, settings)
    if isempty(lines)
        return 0.0
    end
    line = lines[1]
    computed_line = PhotoIonization.computeAmplitudesProperties(line, nm, grid, 100, settings, printout=false)

    return computed_line.crossSection.Coulomb
end
