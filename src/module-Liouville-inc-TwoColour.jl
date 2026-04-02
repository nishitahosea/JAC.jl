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
`Liouville.displayDensityMatrix(stream, liouvilleLevels::Array{Liouville.RamanLevel,1}, densityM::Matrix{ComplexF64})`
    ... Display the current density matrix.
"""
function displayDensityMatrix(stream, liouvilleLevels::Array{Liouville.RamanLevel,1}, densityM::Matrix{ComplexF64})
    println(stream, " ")
    println(stream, "  Selected Liouville levels and current density matrix:")
    println(stream, " ")
    for (idx, liouLevel) in enumerate(liouvilleLevels)
        sa = "       " * string(idx) * ")  "
        sa = sa * string(liouLevel.leadingConfig) * "   "
        sa = sa * string(liouLevel.leadingNotation) * "                          "
        sa = sa[1:70]
        row = densityM[idx, :]
        for z in row
            sa = sa * @sprintf("%8.3f %+8.3fim  ", real(z), imag(z))
        end
        println(sa)
    end
    println(stream, " ")
    return nothing
end
