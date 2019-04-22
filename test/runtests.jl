using Test
using MatrixCalculus

n,m = rand(10:30,2)
n,m = 5,3
A = rand(n,m)
vA = vec(A)
C = comm(n,m,false)

collect(1:m:n*m)
@test C*vec(A) == vec(A')
@test C*vA == vec(A')

A = A*A'
ns = (n+1)n÷2
S = vech(A)
@test length(S) == ns
vech(A,S)
@test S[1] == A[1]
@test S[2] == A[1,2]
@test S[n+1] == A[2,2]
@test S[end] == A[end,end]

E = elim(n)
@test E*vec(A) == vech(A)
