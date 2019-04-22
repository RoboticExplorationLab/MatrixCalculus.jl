module MatrixCalculus
using LinearAlgebra
using SparseArrays
using DocStringExtensions

"$(SIGNATURES) Create a sparse identity matrix"
speye(n::Int) = sparse(Diagonal(1I,n))

"$(SIGNATURES) Create a commutator matrix C_(n,m) such that vec(A') = C*vec(A) where A in an n×m matrix
Will default to a sparse matrix if n*m > 100, or the `sparse` flag is set to `true` "
function comm(n::Int, m::Int, sparse::Bool=false)
    if sparse || n*m > 100
        reshape(kron(speye(n)[:], speye(m)), n*m, n*m)
    else
        reshape(kron(vec(Diagonal(1I,n)), Diagonal(1I,m)), n*m, n*m)
    end
end

"$(SIGNATURES) Returns the half-vectorization of a symmetric matrix as a view of the original matrix.
It is a vector of stacked columns of the lower triangular portion of S"
vech(S) = view(S,trilinds(S))
"$(SIGNATURES) Returns the half-vectorization of a symmetric matrix, given a previous half-vectorization of the same size"
vech(S,S2::SubArray) = view(S,S2.indices[1])

"$(SIGNATURES) Calculate the linear indices of the lower triangular portion of a square matrix"
function trilinds!(inds,A)
    n = size(A,1)
    k = 1
    for i = 1:n
        inds[k:k+n-i] = (i:n) .+ (i-1)n
        k += n-i+1
    end
    inds
end
function trilinds(A)
    n = size(A,1)
    inds = zeros(Int,(n+1)n÷2)
    trilinds!(inds,A)
    return inds
end

"$(SIGNATURES) Calculate the linear indices of the upper triangular portion of a square matrix"
function triuinds!(inds,A)
    n = size(A,1)
    k = 1
    for i = 1:n
        inds[k:k+i-1] = (1:i) .+ (i-1)n
        k += i
    end
    return inds
end
function triuinds(A)
    n = size(A,1)
    inds = zeros(Int,(n+1)n÷2)
    triinds!(inds,A)
    return inds
end

"$(SIGNATURES) Returns the elimination matrix E_n such that vech(S) = E_n*vec(S).
This is useful for taking derivatives with respect to a symmetric matrix"
function elim(n)
    ns = (n+1)n÷2
    E = sparse(1:ns, Array(vech(reshape(1:(n*n),n,n))), ones(Int,ns));
end

"$(SIGNATURES) Returns the duplication matrix D_n such that vec(S) = D_n*vech(S).
This is useful for taking derivatives with respect to a symmetric matrix"
function dupl(n)
    A = ones(Int,n,n)
    A[trilinds(A)] = 1:(n+1)n÷2
    A = Symmetric(A,:L)
    D = sparse(1:(n*n), A[:], ones(Int,n*n))
end

end # module
