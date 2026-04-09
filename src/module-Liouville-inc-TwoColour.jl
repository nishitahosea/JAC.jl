module LiouvilleTwoColour

using ..Basics, ..Defaults, ..Pulse, ..LiouvilleBase
import ..Liouville: AbstractLiouvilleScheme

struct TwoColourScheme <: AbstractLiouvilleScheme
    levelSelection          ::LevelSelection      # Reuse same structure!
    levelNotations          ::Array{String,1}     # Reuse same structure!
    # Add TwoColour-specific parameters if needed
    # e.g., ionizationRate ::Float64
end

function TwoColourScheme()
    TwoColourScheme(LevelSelection(), String[])
end

# Reuse the SAME pattern as StimulatedRamanScheme!
function initializeLevels(scheme::TwoColourScheme, multiplet::Multiplet)
    # Same logic as StimulatedRamanScheme.initializeLevels
    levels = LiouvilleBase.AtomicLevel[]

    for (idx, index) in enumerate(scheme.levelSelection.indices)
        for level in multiplet.levels
            if index == level.index
                leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
                push!(levels, LiouvilleBase.AtomicLevel(leadingConf, scheme.levelNotations[idx], level))
            end
        end
    end

    # Add loss channel for ionization
    if length(scheme.levelNotations) > length(scheme.levelSelection.indices)
        push!(levels, LiouvilleBase.AtomicLevel(Configuration("[He]"), scheme.levelNotations[end], Level()))
    end

    return levels
end

function perform(scheme::TwoColourScheme, computation::Computation; output::Bool=true)
    # Same structure as StimulatedRamanScheme.perform!
    results = Dict{String, Any}()

    # Convert pulses (same as Raman)
    pulses = Pulse.AbstractPulse[]
    for pulse in computation.pulses
        if typeof(pulse) == Pulse.FelPulse
            push!(pulses, Pulse.convertPulse(pulse))
        else
            push!(pulses, pulse)
        end
    end

    # Atomic structure (same as Raman)
    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel,
                                          computation.grid, computation.asfSettings)

    # Initialize levels (same pattern as Raman!)
    levels = initializeLevels(scheme, multiplet)
    densityM = LiouvilleBase.initializeDensityMatrix(levels)  # Same function!

    if computation.settings.printBefore
        LiouvilleBase.displayDensityMatrix(stdout, levels, densityM)  # Same function!
    end

    # HERE is where TwoColour differs - the time evolution!
    # Different Liouville equation, different Hamiltonian, etc.
    println("\n  TwoColour time evolution uses different Liouvillian than Raman.")

    # ... TwoColour-specific time evolution ...

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
    end

    return results
end

end # module
