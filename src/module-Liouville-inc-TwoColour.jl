using ..Basics, ..Defaults, ..Pulse

"""
`struct  Liouville.TwoColourScheme <: Liouville.AbstractLiouvilleScheme`
    ... defines a struct for two-colour XUV+NIR photoionisation.
"""
struct TwoColourScheme <: Liouville.AbstractLiouvilleScheme
    levelSelection          ::LevelSelection
    leadingNotation         ::Array{String,1}
end

# Default constructor
function TwoColourScheme()
    TwoColourScheme(LevelSelection(), String[])
end

# Show function
function Base.show(io::IO, scheme::TwoColourScheme)
    println(io, "TwoColourScheme:")
    println(io, "  levelSelection:    $(scheme.levelSelection)")
    println(io, "  leadingNotation:   $(scheme.leadingNotation)")
end

# Helper function for Gaussian envelope
function gaussianEnvelope(t::Float64, sigma::Float64)
    wa = t^2 / (2*sigma)^2
    return exp(-wa) / (sigma * sqrt(2pi))
end

# Envelope for GaussianSimplified pulse
function envelope(pulse::Pulse.GaussianSimplified, t::Float64)
    sigma = pulse.fwhm / (2 * sqrt(2 * log(2)))  # Convert FWHM to sigma
    wa = (t - pulse.timeDelay)^2 / (2 * sigma^2)
    return exp(-wa) / (sigma * sqrt(2pi))
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
            E_total = E_total + pulse.A0 * envelope(pulse, t) * cos(pulse.omega * (t - pulse.timeDelay) + pulse.cep)
        else
            error("Unknown pulse type: $(typeof(pulse))")
        end
    end
    return E_total
end

# Note: displayDensityMatrix is defined in the Raman file for RamanLevel
# Do NOT redefine it here for TwoColourScheme
