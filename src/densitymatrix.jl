function fermi0(y) 
    
    if y > 0
        return 0.0
    elseif y < 0
        return 1.0
    else
        return 0.5
    end
    
end

function fermi(y, T)
    
    if T == 0 
        return fermi0(y)
    else
        return 0.5 - 0.5*tanh(y/T/2)
    end
    
end

function densitymatrix(ph::PeierlsHamiltonian{E, D, ED, Ti, Ts, Th, Tj}; EF = 0.0, T = 0) where {E, D, ED, Ti, Ts, Th, Tj} 
    
    energies, states = eigen(Hermitian(Matrix(ham(ph))))
    
    N = dimension(ph)
    
    rho = zeros(ComplexF64, size(ph))
    
    for n in eachindex(energies)
        
            rho += fermi(energies[n] - EF, T)*states[:, n]*adjoint(states[:, n])/N
            
    end
    
    return rho
    
end