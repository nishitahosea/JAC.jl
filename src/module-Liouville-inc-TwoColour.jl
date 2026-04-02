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

# DO NOT redefine displayDensityMatrix here — it already exists in the Raman file


"""
`Liouville.perform(scheme::TwoColour, computation::Liouville.Computation; output::Bool=true)`
    ... performs a Liouville time-evolution for two-color XUV+NIR photoionization.
"""
function perform(scheme::TwoColour, computation::Liouville.Computation; output::Bool=true)
    results = Dict{String, Any}()

    println("")
    printstyled("Liouville.perform(): Two-color XUV+NIR computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------ \n", color=:light_green)

    # Step 1: Convert pulses to computational format
    pulses = Pulse.AbstractPulse[]
    for pulse in computation.pulses
        if typeof(pulse) == Pulse.FelPulse
            push!(pulses, Pulse.convertPulse(pulse))
        else
            push!(pulses, pulse)
        end
    end

    # Step 2: Compute atomic structure (SCF)
    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel,
                                          computation.grid, computation.asfSettings)

    # Step 3: Initialize levels, density matrix, and Hamiltonians
    # NOTE: You need to implement these for TwoColour!
    levels = initializeLevels(scheme, multiplet)           # ← need to implement
    densityM = initializeDensityMatrix(levels)             # ← can reuse from Raman
    atomicHM = initializeAtomicHamiltonianMatrix(levels)   # ← can reuse from Raman
    couplingHM = initializeCouplingHamiltonianMatrix(scheme, levels, pulses)  # ← need to implement

    # Step 4: Print initial state
    if computation.settings.printBefore
        displayDensityMatrix(stdout, levels, densityM)
        displayGenericHamiltonian(stdout, levels, atomicHM, couplingHM)  # ← need to implement or adapt
    end

    # Step 5: Time evolution (placeholder for now)
    println("  Time evolution not yet implemented for TwoColour scheme.")
    println("  Initial density matrix has been set up.")

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
        results["atomic_hamiltonian"] = atomicHM
    end

    println("\n> Two-color computation setup complete ...")

    return results
end
