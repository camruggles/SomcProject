

using PyPlot

in = readdlm(ARGS[1])
#n = size(in, 1)
n = parse(Int, ARGS[2])
println("$n")
x = in[1:n, 1]
y = in[1:n, 2]

prime = Array{Float64,1}()
arr = Array{Float64, 1}()

for i in 2:n
	¬ = (y[i] - y[i-1]) / (x[i] - x[i-1])
	push!(prime, ¬)
	push!(arr, x[i])
end

plot(arr, prime)
show()

#=
moved in
got the papers signed
therapy thing
found more datasets
updated my script
derivatives problem
row summation
law of large numbers for one vector

=#
