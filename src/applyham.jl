function applyham!(C::Matrix, A::SparseMatrixCSC, B::Matrix, alpha::Number, beta::Number)
    
    @Threads.threads for j in axes(B, 2)
        @inbounds mul!(view(C, :, j), A, view(B, :, j), alpha, beta)
    end
    
    return C
end

function applyham!(C::Matrix, A::SparseMatrixCSC, B::Matrix)
    
    @Threads.threads for j in axes(B, 2)
        @inbounds mul!(view(C, :, j), A, view(B, :, j))
    end
    
    return C
end

function hermitianpart!(A::Matrix)
    
    n, m = size(A)
    if n != m
        throw(DimensionMismatch("matrix is not square: dimensionas are ($n, $m)"))
    end
    
    for j in 1:m
        
        for i in 1:(j-1)
            @inbounds A[i, j] = conj(A[j, i])
        end
        
        for i in j:n
            @inbounds A[i, j] += conj(A[j, i])
        end
        
    end
    
    return A
    
end

function antihermitianpart!(A::Matrix)
    
    n, m = size(A)
    if n != m
        throw(DimensionMismatch("matrix is not square: dimensionas are ($n, $m)"))
    end
    
    for j in 1:m
        
        for i in 1:(j-1)
            @inbounds A[i, j] = -conj(A[j, i])
        end
        
        for i in j:n
            @inbounds A[i, j] -= conj(A[j, i])
        end
        
    end
    
    return A
    
end