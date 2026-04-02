
using  ..Basics,  ..Defaults, ..Pulse


"""
`struct  Liouville.TwoColour`
    ... defines a struct to comprise the level information for a Liouville time evolution in the two-colour XUV+NIR photoionisation.

    + levelSelection   ::LevelSelection               ... which bound states to include.
    + levelNotations   ::Array{String,1}              ... labels for each level.
"""
struct TwoColour
    levelSelection     ::LevelSelection             # which bound states to include
    leadingNotation    ::Array{String,1}            # labels for each level
    # No frequencies — they come from pulses
    # No includeIonization flag — you always include it if NIR is present
end

"""
`Liouville.TwoColour()`  ... constructor for the default values of Liouville.TwoColour set
"""

# Default constructor
function TwoColour()
    TwoColour(LevelSelection(), String[], false)
end

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################




# `Base.show(io::IO, lioulevel::Liouville.TwoColour)`  ... prepares a proper printout of the variable settings::Liouville.TwoColour.
function Base.show(io::IO, data::Liouville.TwoColour)
    println(io, "TwoColour:")
    println(io, "levelSelection:        $(TwoColour.LevelSelection)  ")
    println(io, "levelNotations:        $(TwoColour.LevelNotations)  ")
end
"""
`Liouville.determineLightField(pulses::Vector{Pulse.AbstractPulse}, t::Float64)`
    ... determines the field amplitude A(t) at the given time t at the (central) position of the atom or
        atomic cloud. All pulses must be of computational type. An fieldAmplitude::Float64 is returned.
"""

function getTotalElectricField(pulses::Vector{Pulse.AbstractPulse}, t::Float64)
    E_total = 0.0

    for pulse in pulses
        if typeof(pulse) == Pulse.GaussianSimplified
            # Check if within pulse duration
            if abs(t - pulse.timeDelay) > 5 * pulse.fwhm
                continue
            end
            # Calculate envelope
            envelope = gaussianEnvelope(t - pulse.timeDelay, pulse.fwhm/2.0)
            # Add contribution with carrier oscillation
            E_total = E_total + pulse.A0 * envelope * cos(pulse.omega * (t - pulse.timeDelay) + pulse.cep)
        else
            error("Unknown pulse type: $(typeof(pulse))")
        end
    end

    return E_total
end

"""
`Liouville.displayDensityMatrix(stream, liouvilleLevels::Array{Liouville.TwoColour,1}, densityM::Matrix{ComplexF64})`
    ... Display the current density matrix together with the ; nothing is returned.
"""
function displayDensityMatrix(stream, liouvilleLevels::Array{Liouville.TwoColour,1}, densityM::Matrix{ComplexF64})

    println(stream, " ")
    println(stream, "  Selected Liouville TwoColour levels and current density matrix:")
    println(stream, " ")
    for  (idx, liouLevel)  in  enumerate(liouvilleLevels)
        sa = "       " * string(idx) * ")  ";   sa = sa[end-10:end]
        sa = sa * string(liouLevel.leadingConfig)    * "   "
        sa = sa * string(liouLevel.leadingNotation)  * "                          "
        sa = sa[1:70]
        row = densityM[idx, :]
        for z in row   sa = sa * @sprintf("%8.3f %+8.3fim  ", real(z), imag(z))    end
        println(sa)
    end
    println(stream, " ")

    return( nothing )
end
