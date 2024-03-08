struct LatticeBasis{E, D, Ts, ED}
    matrix::SMatrix{E, D, Ts, ED}
end

struct HoppingList{Th, Ti}
    rowval::Vector{Ti}
    colval::Vector{Ti}
    nzval::Vector{Th}
end

struct PeierlsHamiltonian{E, D, ED, Ti, Ts, Th, Tj}
    bravais::LatticeBasis{E, D, Ts, ED}
    sites::Vector{SVector{E, Ts}}
    hoppings::HoppingList{Th, Ti}
    separations::Vector{SVector{E, Ts}}
    hamiltonian::SparseMatrixCSC{Th, Ti}
    current::NTuple{E, SparseMatrixCSC{Tj, Ti}}
end

function Base.size(ph::PeierlsHamiltonian)
    
    return (length(ph.sites), length(ph.sites))
    
end

function dimension(ph::PeierlsHamiltonian)
    
    return length(ph.sites)

end

function fractional(x)
    return x - ceil(x-0.5)
end

function pdistance(r, A, Ainv)
    
    u = Ainv*r
    
    return A*fractional.(u)
end

function PeierlsHamiltonian(qh::Quantica.Hamiltonian{T, E, D}) where {T, E, D}
    
    A = SMatrix{E, D}(qh.lattice.bravais.matrix)
    Ainv = pinv(A)
    sites = qh.lattice.unitcell.sites
    
    Is, Js, Vs = findnz(qh[()])
    
    for n in 2:length(qh.harmonics)
        
        In, Jn, Vn = findnz(qh.harmonics[n].h.flat)
        
        append!(Is, In)
        append!(Js, Jn)
        append!(Vs, Vn)
    end
    
    h = sparse(Is, Js, Vs)
    
    rows, cols, hvals = findnz(h)
    
    separations = [pdistance(sites[ij[1]] - sites[ij[2]], A, Ainv) for ij in zip(rows, cols)]

    Jvals = [separations[n]*hvals[n] for n in eachindex(hvals)]

    J = ntuple(i -> sparse(rows, cols, [im*separations[n][i]*hvals[n] for n in eachindex(hvals)]), Val(E))

    return PeierlsHamiltonian(LatticeBasis(A), sites, HoppingList(rows, cols, hvals), separations, h, J)
    
end

function ham(ph::PeierlsHamiltonian, A)
    
    @inbounds for n in eachindex(ph.hamiltonian.nzval)
        
        ph.hamiltonian.nzval[n] = cis(dot(ph.separations[n], A))*ph.hoppings.nzval[n]
        
    end
    
    return ph.hamiltonian
end

function ham(ph::PeierlsHamiltonian)
    
    @inbounds for n in eachindex(ph.hamiltonian.nzval)
        
        ph.hamiltonian.nzval[n] = ph.hoppings.nzval[n]
        
    end
    
    return ph.hamiltonian
end

function current(ph::PeierlsHamiltonian{E, D, ED, Ti, Ts, Th, Tj}, A) where {E, D, ED, Ti, Ts, Th, Tj} 
    
    @inbounds for n in eachindex(ph.hoppings.nzval)
        
        for i in 1:E
            ph.current[i].nzval[n] = im*cis(dot(ph.separations[n], A))*ph.separations[n][i]*ph.hoppings.nzval[n]
        end
        
    end
    
    return ph.current
end

function current(ph::PeierlsHamiltonian{E, D, ED, Ti, Ts, Th, Tj}) where {E, D, ED, Ti, Ts, Th, Tj} 
    
    @inbounds for n in eachindex(ph.hoppings.nzval)
        
        for i in 1:E
            ph.current[i].nzval[n] = im*ph.separations[n][i]*ph.hoppings.nzval[n]
        end
        
    end
    
    return ph.current
end