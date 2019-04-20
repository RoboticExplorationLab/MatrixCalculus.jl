module MatrixCalculus
using LinearAlgebra
using SparseArrays

greet() = print("Hello World!")

speye(n::Int) = sparse(Diagonal(1I,n))
function comm(n::Int, m::Int, sparse::Bool=true)
    if sparse || n*m > 100
        reshape(kron(speye(n)[:], speye(m)), n*m, n*m)
    else
        reshape(kron(vec(Diagonal(1I,n)), Diagonal(1I,m)), n*m, n*m)
    end
end

vech(S) = view(S,triinds(S))
vech(S,S2::SubArray) = view(S,S2.indices[1])

function triinds!(inds,A)
    n = size(A,1)
    k = 1
    for i = 1:n
        inds[k:k+i-1] = (1:i) .+ (i-1)n
        k += i
    end
    return inds
end
function triinds(A)
    n = size(A,1)
    inds = zeros(Int,(n+1)n÷2)
    triinds!(inds,A)
    return inds
end

function elim(n)
    ns = (n+1)n÷2
    E = sparse(1:ns, Array(vech(reshape(1:(n*n),n,n))), ones(Int,ns));
end

function elim2(n)
    ns = (n+1)n÷2
    E = spzeros(Int,ns,n*n)
    r,c = 1,1
    for i = 1:n
        E[r:r+i-1,c:c+i-1] = Diagonal(I,i)
        r += i
        c += n
    end
    return E
end

density(A::SparseMatrixCSC) = nnz(A)/length(A)


end # module
