module LiouvilleTwoColor

using ..Basics, ..Defaults, ..Pulse, ..LiouvilleBase

# Then replace RamanLevel with LiouvilleBase.AtomicLevel
function initializeLevels(scheme::TwoColour, multiplet::Multiplet)
    levels = LiouvilleBase.AtomicLevel[]  # instead of RamanLevel
    # ... rest of code
end

struct TwoColour <: Liouville.AbstractLiouvilleScheme
    levelSelection          ::LevelSelection
    levelNotations          ::Array{String,1}
end

function initializeLevels(scheme::TwoColour, multiplet::Multiplet)
    levels = LiouvilleBase.AtomicLevel[]
    # ... implementation using LiouvilleBase.AtomicLevel
    return levels
end

function perform(scheme::TwoColour, computation::Computation)
    levels = initializeLevels(scheme, multiplet)
    densityM = LiouvilleBase.initializeDensityMatrix(levels)
    LiouvilleBase.displayDensityMatrix(stdout, levels, densityM)
    # ... rest of TwoColor-specific evolution
end

end
