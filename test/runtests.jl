using Test
using BenchmarkTools
using BenchmarkTools

n,m = rand(10:30,2)
A = rand(n,m)
vA = vec(A)
C = comm(n,m,false)
C == comm2(n,m)

collect(1:m:n*m)
@test C*vec(A) == vec(A')
@test C*vA == vec(A')

A = A'A
S = vech(A)
@test length(S) == (m+1)*(m/2)
vech(A,S)
@test S[1] == A[1]
@test S[2] == A[1,2]
@test S[3] == A[2,2]
@test S[end] == A[end,end]


inds = vcat(collect.(([(1:i) .+ (i-1)m for i = 1:m]))...)
(m+1)*m/2
vech2(A)
vech3(A)
S = vech3(A)
@btime vech($A)
@btime vech2($A)
@btime vech3($A)
@btime vech($A,$S)

D = Diagonal(I,10)
sparse(D)
kron(sparsevec(D),sparse(D))

D = Diagonal(I,100)
@btime sparse($D)

@btime spdiagm(0=>ones(100))
@btime sparse(Diagonal(I,100))
D = sparse(Diagonal(I,100))
vec(D)
sparsevec(D1)

D1 = speye(n)
D2 = speye(m)
@btime reshape(kron($D1[:],$D2),n*m,n*m)
@btime reshape(kron(sparsevec($D1),$D2),n*m,n*m)


kron(vec(D1),D2)
kron(sparsevec(D1),D2)
kron(vec(Matrix(1I,n,n)),D2)
sparsevec(D1) == vec(D1)
Matrix(I,n,n)
D1
