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
            if abs(t - pulse.timeDelay) > 5 * pulse.fwhm
                continue
            end
            envelope = gaussianEnvelope(t - pulse.timeDelay, pulse.fwhm/2.0)
            E_total = E_total + pulse.A0 * envelope * cos(pulse.omega * (t - pulse.timeDelay) + pulse.cep)
        else
            error("Unknown pulse type: $(typeof(pulse))")
        end
    end

    return E_total
end

# DO NOT redefine displayDensityMatrix here — it already exists in the Raman file
