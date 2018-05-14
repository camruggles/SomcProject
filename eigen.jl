#z = 3817

using PyPlot
spec = 2828550


#this is a markov chain simulation, the chain makes random jumps to different states and we 
#record the information to see how quickly it becomes accurate to the stationary distribution

function SimulateSOMC(G, station, z, probs)
	counter = 0.0  # this is used to determine how many transitions have been made
	∆ = 0 #this is never used
	a,b = 621, 1478 #arbitrary starting point

	converg = ones(Float64, z)
	instance = zeros(Float64, z)

	iters = [0]
	vals = [0.0]

	indication = true


	for beans in 1:10000000
		 #not sure how to start the chain, doesn't really matter but still
		t = time()
		
		#get graph adj list for state that we're on currenlty
		arr = G[a,b]

		#get a random end state for us to jump to
		i = fRand()
		i = ceil( i/probs[f(a,b,z)] )   #check this line of code
		i = convert(Int64, i)
		#arr = G[a,b]
		c = arr[i]
	
		#update the number of times we've been at a state
		instance[c] += 1
		converg[c] = 1.0 / instance[c]

		#update the convergence distribution each time we get to a state
		indicator = true
		for j in 1:z
			if instance[j] != 0
				#converg[j] = 1.0 / instance[j]
				converg[j] = instance[j] / counter
			else #if it is zero
				indicator = false
			end
		end
		#if there are no zeroes left
		if indicator && indication
			write(STDERR, "All states reached at iteration $beans and at state groupings $a $b $c\n")
			indication = false
		end

		#update the new difference and add it to the x and y axis information
		d = norm(converg-station)
		t = time() - t
		if beans % 10 == 0
			append!(iters, beans)
			append!(vals, d)
			#=
			println("beans: $beans")
			println("\ta: $a, b: $b, c, $c")
			println("\tDiff: $d")
			println("\tTime: $t")
			=#
		end


		a=b
		b=c
		counter += 1

		#this is never used
		if false #d < ∆
			println("below threshold")
			break
		end

	end
	
	println("Iterations until convergence:")
	println("$counter")

	numOnes = 0


	for i in 1:z
		print("My: $(converg[i]) \tStationary: $(station[i])\n")
		if converg[i] == 1
			numOnes += 1
		end
	end
	println("Num ONes: $numOnes")
	println("original norm: $d")

	for i in 1:z
		if instance[i] == 0
			converg[i] = 0
		end
	end

	d = norm(converg-station)
	println("difference in norm: $d")


	return iters, vals
end


#calculates the index of the state pairs
function f(x, y, v)
	(x-1) * v + y
end


#Generates a number between 0 and 1
function fRand()
	s = 2;
	A = abs.(s*rand(1,1) - s/2)
	return A[1]
end

function main()

	A = readdlm("in.txt", Int64)
	z = maximum(A) + 1
	println("Z: $z")
	#reading in probabilities
	probs = Dict{Int64, Float64}()

	P = readdlm("probs.txt", Float64)

	for i in 1:size(P,1 )
		a = convert(Int64, P[i,1]) + 1
		b = convert(Int64, P[i,2]) + 1
		c = P[i,3]
		i1 = f(a,b,z)
		probs[i1] = c

	end

	a = b = c = 0

	#reading in the chain from the text file
	#this is used to store the information about the states for simulation
	G = Dict{Tuple, Array{Int64, 1}}()

	n = size(A,1)
	arg1 = Array{Int64}(n)
	arg2 = Array{Int64}(n)
	arg3 = Array{Float64}(n)

	println("size: $n")

	#intaking input array and then calculating numbers for state pairs
	for i in 1:n
		a = convert(Int64, A[i,1])+1
		b = convert(Int64, A[i,2])+1
		c = convert(Int64, A[i,3])+1
		i1 = f(a,b, z)
		i2 = f(b,c, z)
		arg1[i] = i1
		arg2[i] = i2
		arg3[i] = probs[i1]

		#=
		if i1 == spec
			println("$a, $b, $c")
		end
		=#

		#store extra information in another graph for simulation purposes
		try
			ar = G[a,b]
			append!(ar, [c])
			ar=0
		catch
			G[a,b] = [c]
		end
	end


	#constructing the matrix
	println("Constructing matrix")
	@time M = sparse(arg2, arg1, arg3)

	m = minimum(M)
	println("min: $m")
	#quit()

	#free variables
	arg2 = arg1 = arg3 = P = A = 0
	gc()

	#getting information about the size of the eigenvectors
	n = size(M,1)
	tmp = size(M,2)
	println("Sizes: $n, $tmp")



	#confirming the distribution is correct
#=
	ei = ones(1,n)
	tmp = M[1:n, spec]
	cameron = ei * tmp
	println("$cameron")
=#

#=
	cam = probs[spec]
	println("probs count: $cam")


	#742, 153, 2570
	i1 = f(742, 153, 3817)
	i2 = f(153, 2570, 3817)
	i = M[i2, i1]
	println("742, 153, 2570 index")
	println("$i")
	quit()
	=#

	#confirming stochasticity of the matrix
	#is not row stochastic
	#confirmed is column stochastic


	#=
	stoch = ei*M
	for i in stoch
		if abs(i-1.0) < 0.01 || abs(i-0.0) < 0.01
			println("Good: $i")
		else
			println("Problem: $i")
		end

	end
	println("end stochasticity testing")

	quit()
	=#


#=
	#Resume normal funcitonality
	#calculated, isolate, and print the eigenvalues and eigenvector information
	println("Eigenvector calculation")
	@time eigenvalues,eigenvectors = eigs(M, nev=5)


	eigenSize = size(eigenvectors, 1)
	c = eigenvectors[1:eigenSize, 1]
	str = typeof(c)
	println("Type: $str")


	println("Eigenvalues")
	println("$eigenvalues")

	println("Eigenvector Size")
	println("$eigenSize")


	#computing the sum using the inner product
	#again not useful
	ei = ones(1, eigenSize)
	inner = ei*c

	println("inner product : $inner")

	#doing the same thing but in a more primitive way
	s = 0.0
	for i in c
		s += i
	end

	println("sum: $s")

	c = c/s


	#calculating the one state stationary distribution, instead of over pairs of states
	station = zeros(Float64, z, 1)
	for i in 1:z
		for j in i:z:eigenSize
			station[i] += c[j]
		end
	end

	#sum of the stationary distribution
	#not useful except for error checking
	#should be one
	ei = ones(1,z)
	inner = ei*station
	println("Stationary sum")
	println("$inner")


	#this looks at the initial difference between the stationary distro and a vector of ones
	ei = ones(z,1)
	@time d = norm(ei-station)
	println("d: $d")
	#quit()
=#
	#pretty straight forward description
	station = readdlm("stationaryDist.out")
	write(STDERR, "beginning simulation\n")

	#begin simulation, and then plot the convergence rate and print axis information to file
	@time iters, vals = SimulateSOMC(G, station, z, probs)
	help = [iters vals]
	writedlm("infoOut.txt", help)
	plot(iters, vals)
	show()

end

write(STDERR, "begin\n")
main()