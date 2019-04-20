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

vech(S) = view(S,UpperTriangular(ones(Bool,size(S))))
vech(S,S2::SubArray) = view(S,S2.indices[1])

function vech2(S)
    n = size(S,1)
    inds = vcat(collect.(([(1:i) .+ (i-1)n for i = 1:n]))...)
    view(S,inds)
end

function vech3(S)
    n = size(S,1)
    inds = zeros(Int,(n+1)n÷2)
    k = 0
    for i = 1:n
        inds[k:k+i] = (1:i) .+ (i-1)n
        k += i
    end
    view(S,inds)
end

function comm2(n,m,sparse::Bool=false)
    if sparse || n*m > 100
        C = spzeros(n*m,n*m)
    else
        C = zeros(Int,n*m,n*m)
    end
    for j = 1:m
        inds = (1:m:n*m) .+ (j-1)
        for i = 1:n
            C[inds[i],i+n*(j-1)] = 1
        end
    end
    return C
end


end # module
