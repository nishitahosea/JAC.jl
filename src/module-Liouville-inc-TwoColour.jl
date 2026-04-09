# module-Liouville-inc-TwoColour.jl
# TwoColour scheme for XUV+NIR photoionization - following StimulatedRamanScheme pattern

using Printf
using ..Basics, ..Defaults, ..Pulse, ..LiouvilleBase

# Import the abstract type from the parent Liouville module
import ..Liouville: AbstractLiouvilleScheme

"""
`struct TwoColourScheme <: AbstractLiouvilleScheme`
    ... defines a scheme for two-color (XUV+NIR) photoionization computations.
    Follows the same pattern as StimulatedRamanScheme.
"""
struct TwoColourScheme <: AbstractLiouvilleScheme
    levelSelection          ::LevelSelection      # Which levels to include
    levelNotations          ::Array{String,1}     # Names for each level
    # TwoColour-specific parameters can be added here
    # e.g., ionizationRate   ::Float64
end

"""
`TwoColourScheme()` - Default constructor
"""
function TwoColourScheme()
    TwoColourScheme(LevelSelection(), String[])
end

"""
`Base.show(io::IO, scheme::TwoColourScheme)` - Print scheme information
"""
function Base.show(io::IO, scheme::TwoColourScheme)
    println(io, "TwoColourScheme:")
    println(io, "  levelSelection:    $(scheme.levelSelection)")
    println(io, "  levelNotations:    $(scheme.levelNotations)")
end

"""
`initializeLevels(scheme::TwoColourScheme, multiplet::Multiplet)`
    - Initialize atomic levels following the same pattern as StimulatedRamanScheme
"""
function initializeLevels(scheme::TwoColourScheme, multiplet::Multiplet)
    levels = LiouvilleBase.AtomicLevel[]

    # Add selected levels
    for (idx, index) in enumerate(scheme.levelSelection.indices)
        for level in multiplet.levels
            if index == level.index
                leadingConf = Basics.extractConfiguration(Basics.LeadingConfiguration(), level)
                push!(levels, LiouvilleBase.AtomicLevel(leadingConf, scheme.levelNotations[idx], level))
            end
        end
    end

    # Add loss channel (ionization continuum)
    if length(scheme.levelNotations) > length(scheme.levelSelection.indices)
        push!(levels, LiouvilleBase.AtomicLevel(Configuration("[He]"), scheme.levelNotations[end], Level()))
    end

    return levels
end

"""
`perform(scheme::TwoColourScheme, computation::Computation; output::Bool=true)`
    - Perform two-color XUV+NIR photoionization computation
"""
function perform(scheme::TwoColourScheme, computation::Computation; output::Bool=true)
    results = Dict{String, Any}()

    println("")
    printstyled("Liouville.perform(): TwoColour XUV+NIR computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------ \n", color=:light_green)

    # Convert pulses (same as StimulatedRamanScheme)
    pulses = Pulse.AbstractPulse[]
    for pulse in computation.pulses
        if typeof(pulse) == Pulse.FelPulse
            push!(pulses, Pulse.convertPulse(pulse))
        else
            push!(pulses, pulse)
        end
    end

    # Atomic structure (same as StimulatedRamanScheme)
    multiplet = SelfConsistent.performSCF(computation.refConfigs, computation.nuclearModel,
                                          computation.grid, computation.asfSettings)

    # Initialize levels (same pattern as StimulatedRamanScheme)
    levels = initializeLevels(scheme, multiplet)
    densityM = LiouvilleBase.initializeDensityMatrix(levels)

    # Print initial state if requested (same as StimulatedRamanScheme)
    if computation.settings.printBefore
        LiouvilleBase.displayDensityMatrix(stdout, levels, densityM)
    end

    println("\n  Time evolution not yet implemented for TwoColourScheme.")
    println("  (Will use different Liouvillian than StimulatedRamanScheme)")

    if output
        results["levels"] = levels
        results["initial_density"] = densityM
    end

    return results
end
