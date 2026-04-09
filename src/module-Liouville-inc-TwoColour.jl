# module-Liouville-inc-TwoColour.jl
# NO module declaration at the top!

using Printf
using ..Basics, ..Defaults, ..Pulse, ..LiouvilleBase
import ..Liouville: AbstractLiouvilleScheme

struct TwoColourScheme <: AbstractLiouvilleScheme
    levelSelection          ::LevelSelection
    levelNotations          ::Array{String,1}
end

function TwoColourScheme()
    TwoColourScheme(LevelSelection(), String[])
end

function Base.show(io::IO, scheme::TwoColourScheme)
    println(io, "TwoColourScheme:")
    println(io, "  levelSelection:    $(scheme.levelSelection)")
    println(io, "  levelNotations:    $(scheme.levelNotations)")
end

function initializeLevels(scheme::TwoColourScheme, multiplet::Multiplet)
    levels = LiouvilleBase.AtomicLevel[]

    for (idx, index) in enumerate(scheme.levelSelection.indices)
        for level in multiplet.levels
            if index == level.index
                leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
                push!(levels, LiouvilleBase.AtomicLevel(leadingConf, scheme.levelNotations[idx], level))
            end
        end
    end

    if length(scheme.levelNotations) > length(scheme.levelSelection.indices)
        push!(levels, LiouvilleBase.AtomicLevel(Configuration("[He]"), scheme.levelNotations[end], Level()))
    end

    return levels
end

function perform(scheme::TwoColourScheme, computation::Computation; output::Bool=true)
    results = Dict{String, Any}()

    println("")
    printstyled("Liouville.perform(): TwoColour XUV+NIR computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------ \n", color=:light_green)

    pulses = Pulse.AbstractPulse[]
    for pulse in computation.pulses
        if typeof(pulse) == Pulse.FelPulse
            push!(pulses, Pulse.convertPulse(pulse))
        else
            push!(pulses, pulse)
        end
    end

    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel,
                                          computation.grid, computation.asfSettings)

    levels = initializeLevels(scheme, multiplet)
    densityM = LiouvilleBase.initializeDensityMatrix(levels)

    if computation.settings.printBefore
        LiouvilleBase.displayDensityMatrix(stdout, levels, densityM)
    end

    println("\n  TwoColour time evolution not yet implemented.")

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
    end

    return results
end

# NO 'end' for a module here!
