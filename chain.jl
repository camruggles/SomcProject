
function fRand()
	s = 2;
	A = abs.(s*rand(1,1) - s/2)
	return A[1]
end




G = Dict{Int64, Array{Int64, 1}}()
inst = Dict{Int64, Int64}()
inst[1] = inst[2] = inst[3] = 0
c = 0.0

G[1] = [1, 1, 2, 3]
G[2] = [1, 3]
G[3] = [1, 2, 3, 3]

state = 1
inst[state] += 1
station = zeros(Float64, 3)
station[state] = 1.0/inst[state]
eigenvector = [0.4; 0.2; 0.4]

n = 10000
diffs = zeros(Float64, n, 1)
indices = zeros(Int64, n, 1)

for times in 1:n
	arr = G[state]


	#println("from state: $state")

	i = fRand()
	#println("i: $i")
	#i = ceil( i/probs[f(a,b,z)] )
	prob = 1.0 / size(arr, 1)
	#println("prob: $prob")
	i = ceil( i/prob )
	i = convert(Int64, i)
	#println("i: $i")
	#arr = G[a,b]
	state = arr[i]


	#println("to state : $state")
	c += 1
	inst[state] += 1
	station[state] = 1/inst[state]

	#v1, v2, v3 = inst[1]/c, inst[2]/c, inst[3]/c

	for j in 1:length(station)
#		station[j] = inst[j] / c
	end

	println("$station")
	#print("\n\n")



	d = norm(station-eigenvector)
	indices[times] = times
	diffs[times] = d
	


end
println("$station")

using PyPlot
plot(indices, diffs)
show()


