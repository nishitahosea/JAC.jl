
using  ..Basics,  ..Defaults, ..Pulse


"""
`struct  Liouville.TwoColourScheme`
    ... defines a struct to comprise the level information for a Liouville time evolution in the two-colour XUV+NIR photoionisation.

    + levelSelection   ::LevelSelection               ... which bound states to include.
    + levelNotations   ::Array{String,1}              ... labels for each level.
"""
struct TwoColourScheme
    levelSelection     ::LevelSelection             # which bound states to include
    leadingNotation    ::Array{String,1}            # labels for each level
    # No frequencies — they come from pulses
    # No includeIonization flag — you always include it if NIR is present
end

# Default constructor
function TwoColourScheme()
    TwoColourScheme(LevelSelection(), String[], false)
end

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`Liouville.TwoColourScheme()`  ... constructor for the default values of Liouville.TwoColour set
"""
function TwoColourScheme()
    TwoColourScheme( LevelSelection(), Level() )
end


# `Base.show(io::IO, lioulevel::Liouville.TwoColourScheme)`  ... prepares a proper printout of the variable settings::Liouville.TwoColourScheme.
function Base.show(io::IO, data::Liouville.RamanLevel)
    println(io, "TwoColourScheme:")
    println(io, "levelSelection:        $(scheme.LevelSelection)  ")
    println(io, "levelNotations:        $(scheme.LevelNotations)  ")
end
