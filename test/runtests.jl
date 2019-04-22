using MatrixCalculus
using Test
using SparseArrays

# Commutator matrix
n,m = rand(11:30,2)
A = rand(n,m)
vA = vec(A)
C = comm(n,m,false)
@test issparse(C.parent)
@test !(C isa Array)
@test C*vec(A) == vec(A')
@test C*vA == vec(A')

C2 = comm(4,4,false)
@test C2 isa Array
C2 = comm(4,4,true)
@test !(C2 isa Array)
@test issparse(C2.parent)

# vech
A = A*A'
ns = (n+1)n÷2
S = vech(A)
@test length(S) == ns
vech(A,S)
@test S[1] == A[1]
@test S[2] == A[1,2]
@test S[n+1] == A[2,2]
@test S[end] == A[end,end]

# Triangular indices
B = rand(n,n)
linds = MatrixCalculus.trilinds(B)
MatrixCalculus.trilinds!(linds,B)
@test linds[1] == 1
@test linds[2] == 2
@test linds[n] == n
@test linds[n+1] == n+2
@test linds[end] == n^2
@test length(linds) == ns

uinds = MatrixCalculus.triuinds(B)
MatrixCalculus.triuinds!(uinds,B)
@test uinds[1] == 1
@test uinds[2] == n+1
@test uinds[3] == n+2
@test uinds[end] == n^2
@test length(uinds) == ns

# Elim matrix
E = elim(n)
@test E*vec(A) == vech(A)

# Dupl matrix
D = dupl(n)
@test vec(A) == D*vech(A)
