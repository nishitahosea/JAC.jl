using ..Basics, ..Defaults, ..Pulse

# Helper function for Gaussian envelope
function gaussianEnvelope(t::Float64, sigma::Float64)
    wa = t^2 / (2*sigma)^2
    return exp(-wa) / (sigma * sqrt(2pi))
end

"""
`struct  Liouville.TwoColour <: Liouville.AbstractLiouvilleScheme`
    ... defines a struct for two-colour XUV+NIR photoionisation.
"""
struct TwoColour <: Liouville.AbstractLiouvilleScheme
    levelSelection          ::LevelSelection
    levelNotations          ::Array{String,1}
end

# Default constructor
function TwoColour()
    TwoColour(LevelSelection(), String[])
end

# Show function
function Base.show(io::IO, scheme::TwoColour)
    println(io, "TwoColour:")
    println(io, "  levelSelection:    $(scheme.levelSelection)")
    println(io, "  levelNotations:    $(scheme.levelNotations)")
end

"""
`Liouville.initializeLevels(scheme::TwoColour, multiplet::Multiplet)`
    ... selects levels based on scheme.levelSelection.
"""
function initializeLevels(scheme::TwoColour, multiplet::Multiplet)
    liouvilleLevels = Liouville.RamanLevel[]
    noLevels = length(scheme.levelSelection.indices)
    foundIndices = falses(noLevels)

    if length(scheme.levelNotations) != noLevels + 1
        error("Expect $(noLevels + 1) strings for level notations (including loss channel), got $(length(scheme.levelNotations))")
    end

    for (idx, index) in enumerate(scheme.levelSelection.indices)
        for level in multiplet.levels
            if index == level.index
                foundIndices[idx] = true
                leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
                liouvLevel = Liouville.RamanLevel(leadingConf, scheme.levelNotations[idx], level)
                push!(liouvilleLevels, liouvLevel)
            end
        end
    end

    if !all(foundIndices)
        error("Not all level indices found in multiplet; foundIndices = $foundIndices")
    end

    # Add loss channel
    push!(liouvilleLevels, Liouville.RamanLevel(Configuration("[He]"), scheme.levelNotations[end], Level()))

    return liouvilleLevels
end

"""
`Liouville.initializeDensityMatrix(liouvilleLevels::Array{Liouville.RamanLevel,1})`
    ... initializes density matrix with ground state populated.
"""
function initializeDensityMatrix(liouvilleLevels::Array{Liouville.RamanLevel,1})
    noLevels = length(liouvilleLevels)
    densityM = zeros(ComplexF64, noLevels, noLevels)
    energies = [lioulevel.level.energy for lioulevel in liouvilleLevels]

    lowestEn = minimum(energies)
    idx = findfirst(==(lowestEn), energies)
    densityM[idx, idx] = 1.0

    return densityM
end

"""
`Liouville.initializeAtomicHamiltonianMatrix(scheme::TwoColour, liouvilleLevels::Array{Liouville.RamanLevel,1})`
    ... initializes the atomic Hamiltonian matrix.
"""
function initializeAtomicHamiltonianMatrix(scheme::TwoColour, liouvilleLevels::Array{Liouville.RamanLevel,1})
    noLevels = length(liouvilleLevels)
    energies = Float64[]
    hamiltonian = zeros(ComplexF64, noLevels, noLevels)

    for lioulevel in liouvilleLevels
        push!(energies, lioulevel.level.energy)
    end
    lowestEn = minimum(energies)

    for n in 1:noLevels
        hamiltonian[n,n] = energies[n] - lowestEn
    end

    return hamiltonian
end

"""
`Liouville.initializeCouplingHamiltonianMatrix(scheme::TwoColour, levels::Array{Liouville.RamanLevel,1}, pulses)`
    ... placeholder for V(t) coupling.
"""
function initializeCouplingHamiltonianMatrix(scheme::TwoColour, levels::Array{Liouville.RamanLevel,1}, pulses)
    noLevels = length(levels)
    couplingHM = Matrix{Function}(undef, noLevels, noLevels)
    for i in 1:noLevels, j in 1:noLevels
        couplingHM[i,j] = t -> 0.0 + 0.0im
    end
    return couplingHM
end

"""
`Liouville.getTotalElectricField(pulses::Vector{Pulse.AbstractPulse}, t::Float64)`
    ... determines the total electric field E(t) at time t.
"""
function getTotalElectricField(pulses::Vector{Pulse.AbstractPulse}, t::Float64)
    E_total = 0.0
    for pulse in pulses
        if typeof(pulse) == Pulse.GaussianSimplified
            if abs(t - pulse.timeDelay) > 5 * pulse.fwhm
                continue
            end
            envelope = gaussianEnvelope(t - pulse.timeDelay, pulse.fwhm/2.0)
            E_total = E_total + pulse.A0 * envelope * cos(pulse.omega * (t - pulse.timeDelay) + pulse.cep)
        else
            error("Unknown pulse type: $(typeof(pulse))")
        end
    end
    return E_total
end

"""
`Liouville.perform(scheme::TwoColour, computation::Liouville.Computation; output::Bool=true)`
    ... main driver for two-color computation.
"""
function perform(scheme::TwoColour, computation::Liouville.Computation; output::Bool=true)
    results = Dict{String, Any}()

    println("")
    printstyled("Liouville.perform(): Two-color XUV+NIR computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------ \n", color=:light_green)

    # Step 1: Convert pulses
    pulses = Pulse.AbstractPulse[]
    for pulse in computation.pulses
        if typeof(pulse) == Pulse.FelPulse
            compPulse = Pulse.convertPulse(pulse)
            push!(pulses, compPulse)
        else
            push!(pulses, pulse)
        end
    end

    # Step 2: Atomic structure
    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel,
                                          computation.grid, computation.asfSettings)

    # Step 3: Initialize
    levels = initializeLevels(scheme, multiplet)
    densityM = initializeDensityMatrix(levels)
    atomicHM = initializeAtomicHamiltonianMatrix(scheme, levels)
    couplingHM = initializeCouplingHamiltonianMatrix(scheme, levels, pulses)

    # Step 4: Print
    if computation.settings.printBefore
        displayDensityMatrix(stdout, levels, densityM)
        println("\n  Atomic Hamiltonian (diagonal):")
        for i in 1:size(atomicHM,1)
            println("    H[$i,$i] = $(real(atomicHM[i,i])) a.u.")
        end
    end

    println("\n  Time evolution not yet implemented for TwoColour scheme.")
    println("  Initial density matrix has been set up.")

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
        results["atomic_hamiltonian"] = atomicHM
    end

    println("\n> Two-color computation setup complete ...")

    return results
end
