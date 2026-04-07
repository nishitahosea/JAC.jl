module LiouvilleStimulatedRaman

using ..Basics, ..Defaults, ..Pulse, ..LiouvilleBase

struct StimulatedRamanScheme <: Liouville.AbstractLiouvilleScheme
    levelSelection          ::LevelSelection
    levelNotations          ::Array{String,1}
    gammaR                  ::Float64
    gammaA                  ::Float64
    calcPopulations         ::Bool
end

function initializeLevels(scheme::StimulatedRamanScheme, multiplet::Multiplet)
    levels = LiouvilleBase.AtomicLevel[]
    # ... implementation using LiouvilleBase.AtomicLevel
    return levels
end

function perform(scheme::StimulatedRamanScheme, computation::Computation)
    levels = initializeLevels(scheme, multiplet)
    densityM = LiouvilleBase.initializeDensityMatrix(levels)
    LiouvilleBase.displayDensityMatrix(stdout, levels, densityM)
    # ... rest of Raman-specific evolution
end

end
