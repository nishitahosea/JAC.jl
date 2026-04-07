using ..Basics, ..Defaults, ..Pulse

"""
`struct  Liouville.TwoColour <: Liouville.AbstractLiouvilleScheme`
    ... defines a struct for two-colour XUV+NIR photoionisation.
"""
struct TwoColour <: Liouville.AbstractLiouvilleScheme
    levelSelection          ::LevelSelection
    leadingNotation         ::Array{String,1}
end

# Default constructor
function TwoColour()
    TwoColour(LevelSelection(), String[])
end

# Show function
function Base.show(io::IO, scheme::TwoColour)
    println(io, "TwoColour:")
    println(io, "  levelSelection:    $(scheme.levelSelection)")
    println(io, "  leadingNotation:   $(scheme.leadingNotation)")
end

# Helper function for Gaussian envelope
function gaussianEnvelope(t::Float64, sigma::Float64)
    wa = t^2 / (2*sigma)^2
    return exp(-wa) / (sigma * sqrt(2pi))
end

# Envelope for GaussianSimplified pulse
function envelope(pulse::Pulse.GaussianSimplified, t::Float64)
    sigma = pulse.fwhm / (2 * sqrt(2 * log(2)))
    wa = (t - pulse.timeDelay)^2 / (2 * sigma^2)
    return exp(-wa) / (sigma * sqrt(2pi))
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
            E_total = E_total + pulse.A0 * envelope(pulse, t) * cos(pulse.omega * (t - pulse.timeDelay) + pulse.cep)
        else
            error("Unknown pulse type: $(typeof(pulse))")
        end
    end
    return E_total
end

"""
`Liouville.perform(scheme::TwoColour, computation::Liouville.Computation; output::Bool=true)`
    ... performs a Liouville time-evolution for two-color XUV+NIR photoionization.
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

    # Step 3: Initialize levels
    levels = Liouville.RamanLevel[]
    for (idx, notation) in enumerate(scheme.leadingNotation)
        if idx <= length(multiplet.levels)
            level = multiplet.levels[idx]
            leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
            push!(levels, Liouville.RamanLevel(leadingConf, notation, level))
        else
            push!(levels, Liouville.RamanLevel(Configuration("[He]"), notation, Level()))
        end
    end

    # Initialize density matrix
    densityM = zeros(ComplexF64, length(levels), length(levels))
    minEnergy = Inf
    minIdx = 1
    for (idx, lvl) in enumerate(levels)
        if lvl.level.energy < minEnergy && lvl.level.energy != 0.0
            minEnergy = lvl.level.energy
            minIdx = idx
        end
    end
    densityM[minIdx, minIdx] = 1.0

    # Print initial state
    if computation.settings.printBefore
        println("\n  Selected Liouville levels and current density matrix:")
        for (idx, lvl) in enumerate(levels)
            pop = real(densityM[idx, idx])
            println("    Level $idx: $(lvl.leadingNotation) → Population = $(@sprintf("%.6f", pop))")
        end
    end

    println("\n  Time evolution not yet implemented for TwoColour scheme.")
    println("  Initial density matrix has been set up.")

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
    end

    println("\n> Two-color computation setup complete ...")

    return results
end
