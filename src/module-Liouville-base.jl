module LiouvilleBase

using ..Basics, ..Defaults, ..Pulse

# Common level container (renamed from RamanLevel)
struct AtomicLevel
    leadingConfig      ::Configuration 
    leadingNotation    ::String  
    level              ::Level  
end

# Common function to initialize density matrix (ground state populated)
function initializeDensityMatrix(levels::Vector{AtomicLevel})
    noLevels = length(levels)
    densityM = zeros(ComplexF64, noLevels, noLevels)
    energies = [lvl.level.energy for lvl in levels]
    
    # Find ground state (lowest energy, excluding loss channel placeholder)
    valid_energies = [e for e in energies if e != 0.0]
    if isempty(valid_energies)
        idx = 1  # fallback
    else
        lowestEn = minimum(valid_energies)
        idx = findfirst(==(lowestEn), energies)
    end
    densityM[idx, idx] = 1.0
    
    return densityM
end

# Common function to display density matrix
function displayDensityMatrix(stream, levels::Vector{AtomicLevel}, densityM::Matrix{ComplexF64})
    println(stream, " ")
    println(stream, "  Selected levels and current density matrix:")
    println(stream, " ")
    for (idx, lvl) in enumerate(levels)
        pop = real(densityM[idx, idx])
        println(stream, "    Level $idx: $(lvl.leadingNotation) → Population = $(@sprintf("%.6f", pop))")
    end
    println(stream, " ")
    return nothing
end

# Common function to initialize atomic Hamiltonian
function initializeAtomicHamiltonian(levels::Vector{AtomicLevel})
    noLevels = length(levels)
    H = zeros(ComplexF64, noLevels, noLevels)
    energies = [lvl.level.energy for lvl in levels]
    minEnergy = minimum([e for e in energies if e != 0.0])
    
    for i in 1:noLevels
        H[i,i] = energies[i] - minEnergy
    end
    return H
end

end # module