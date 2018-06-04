
networkName = ARGS[1]
A = readdlm("$networkName-In.txt", Int64)

for i in eachindex(A)
  A[i] += 1
end

D = Dict{Int64, Int64}()
S = IntSet(A)
c = 0



for i in S
  println("$i : $c")
  D[i] = c
  c += 1

end


s1 = size(A, 1)
s2 = size(A, 2)

for i in 1:s1
   A[i, 1] = D[A[i, 1]]
   A[i, 2] = D[A[i, 2]]
   A[i, 3] = D[A[i, 3]]
   
end

writedlm("UpdatedInts.txt", A)

A = readdlm("$networkName-Probs.txt", Float64)
n = size(A, 1)
for i in 1:n

   A[i, 1] += 1
   A[i, 2] += 1

   A[i, 1] = D[A[i, 1]]
   A[i, 2] = D[A[i, 2]]
end

writedlm("ProbsUpdated.txt", A)

#=
propose a new approach to solve an existing problem
test that approach with experiments
because most algorithms won't just involve quick analysis with asymptotic complexity
=#
