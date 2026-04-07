module LiouvilleStimulatedRaman

using Printf
using ..Basics, ..Defaults, ..Pulse, ..LiouvilleBase

# Import the abstract type from the parent Liouville module
import ..Liouville: AbstractLiouvilleScheme

struct StimulatedRamanScheme <: AbstractLiouvilleScheme
    levelSelection          ::LevelSelection
    levelNotations          ::Array{String,1}
    gammaR                  ::Float64
    gammaA                  ::Float64
    calcPopulations         ::Bool
end

function StimulatedRamanScheme()
    StimulatedRamanScheme(LevelSelection(), String[], 0.0, 0.0, false)
end

function Base.show(io::IO, scheme::StimulatedRamanScheme)
    println(io, "StimulatedRamanScheme:")
    println(io, "  levelSelection:    $(scheme.levelSelection)")
    println(io, "  levelNotations:    $(scheme.levelNotations)")
    println(io, "  gammaR:            $(scheme.gammaR)")
    println(io, "  gammaA:            $(scheme.gammaA)")
    println(io, "  calcPopulations:   $(scheme.calcPopulations)")
end

# Initialize levels using LiouvilleBase.AtomicLevel
function initializeLevels(scheme::StimulatedRamanScheme, multiplet::Multiplet)
    levels = LiouvilleBase.AtomicLevel[]
    noLevels = length(scheme.levelSelection.indices)

    for (idx, index) in enumerate(scheme.levelSelection.indices)
        for level in multiplet.levels
            if index == level.index
                leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
                push!(levels, LiouvilleBase.AtomicLevel(leadingConf, scheme.levelNotations[idx], level))
            end
        end
    end

    # Add loss channel
    push!(levels, LiouvilleBase.AtomicLevel(Configuration("[He]"), scheme.levelNotations[end], Level()))

    return levels
end

# Perform function for Raman
function perform(scheme::StimulatedRamanScheme, computation::Computation; output::Bool=true)
    results = Dict{String, Any}()

    println("")
    printstyled("Liouville.perform(): Stimulated Raman computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------ \n", color=:light_green)

    # Convert pulses
    pulses = Pulse.AbstractPulse[]
    for pulse in computation.pulses
        if typeof(pulse) == Pulse.FelPulse
            push!(pulses, Pulse.convertPulse(pulse))
        else
            push!(pulses, pulse)
        end
    end

    # Atomic structure
    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel,
                                          computation.grid, computation.asfSettings)

    # Initialize levels
    levels = initializeLevels(scheme, multiplet)
    densityM = LiouvilleBase.initializeDensityMatrix(levels)

    if computation.settings.printBefore
        LiouvilleBase.displayDensityMatrix(stdout, levels, densityM)
    end

    println("\n  Time evolution not yet fully implemented for StimulatedRaman scheme.")

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
    end

    return results
end

end # module
