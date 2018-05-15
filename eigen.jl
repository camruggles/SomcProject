#TODO
#=
#Massively Refactor Code

Test the norms between eigenvectors of successive calculations, and for the stationary distributions as well

=#


using PyPlot

networkName = "NAN"
stateA = 621
stateB = 1478


function extractStationaryDist(maxState, c, eigenSize)
	station = zeros(Float64, maxState, 1)
	for i in 1:maxState
		for j in i:maxState:eigenSize
			station[i] += c[j]
		end
	end

	#sum of the stationary distribution
	#not useful except for error checking
	#should be one
	ei = ones(1, maxState)
	inner = ei*station
	println("Stationary sum")
	println("$inner")
	if abs(inner[1] - 1.0) > 0.01
		write(STDERR, "Error : Stationary Distribution does not add to one")
		quit()
	end

	return station
end

function extractEigenvector(M)
	#calculated, isolate, and print the eigenvalues and eigenvector information
	println("Eigenvector calculation")
	@time eigenvalues,eigenvectors = eigs(M, nev=5)

	eigen_id = -1
	for i in 1:5
		println("Eigenvalue $i $(eigenvalues[i])")
		if abs(eigenvalues[i] - 1.0+0.0im) < 0.001
			eigen_id = i
			break
		end
	end

	if eigen_id == -1
		write(STDERR, "Error, no eigenvector with eigenvalue of one")
		quit()
	end




	eigenSize = size(eigenvectors, 1)
	c = eigenvectors[1:eigenSize, eigen_id]

	println("Eigenvalues")
	println("$eigenvalues")

	println("Eigenvector Size")
	println("$eigenSize")


	c = convert(Array{Float64, 1}, c)


	return c, eigenSize
end

function extractNormalizedEigenvector(M)

	c, eigenSize = extractEigenvector(M)

	ei = ones(1, eigenSize)
	inner = ei*c

	println("inner product : $inner")
	s = inner[1]
	c = c/s

	return c, eigenSize
end

function getMatrix(inputFile, probsFile)
	A = readdlm(inputFile, Int64)
	maxState = maximum(A) + 1
	println("Max State: $maxState")
	#reading in probabilities
	probs = Dict{Int64, Float64}()

	P = readdlm(probsFile, Float64)
	for i in 1:size(P,1 )
		a = convert(Int64, P[i,1]) + 1
		b = convert(Int64, P[i,2]) + 1
		c = P[i,3]
		i1 = f(a,b,maxState)
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

	#i j k
	#i1 = i j
	#i2 = j k
	#probs is ij but should be jk
	#from jk it has n destinations i with P(1/n)
	#instead, we do ij to destination k with P(1/ij)
	#this works as if it were also i <- jk as well as ij->k
	#since the state ordering is the relevant part and it's super symmetric


	psize = size(A, 1)

	inl = psize * 4 / 10

	inl = round(Int64, inl)
	global stateA = convert(Int64, A[inl,1]) + 1
	global stateB = convert(Int64, A[inl,2]) + 1



	for i in 1:n
		#intaking the 3 variables
		a = convert(Int64, A[i,1])+1
		b = convert(Int64, A[i,2])+1
		c = convert(Int64, A[i,3])+1

		#converting the indices to states
		i1 = f(a,b, maxState)
		i2 = f(b,c, maxState)

		#recording the information in arrays for quicker Spare matrix construction
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

	#free variables
	arg2 = arg1 = arg3 = P = A = 0
	gc()

	return G, M, maxState, probs
end	

function testColumnStochastic(M, n)

	
	ei = ones(1,n)
	stoch = ei*M
	stochastic = true
	for i in stoch

		if abs(i-1.0) < 0.01 || abs(i-0.0) < 0.01
			#println("Good: $i")
		else
			#println("Problem: $i")
			stochastic = false
		end

	end

	println("end stochasticity testing")
	println("Column Stochastic: $stochastic")
	return stochastic
end


#this is a markov chain simulation, the chain makes random jumps to different states and we 
#record the information to see how quickly it becomes accurate to the stationary distribution

function SimulateSOMC(G, station, maxState, probs, oneVec = false)
	counter = 0.0  # this is used to determine how many transitions have been made
	∆ = 0 #this is never used
	a,b = stateA, stateB #arbitrary starting point

	converg = 0
	if oneVec
		converg = ones(Float64, maxState)
	else
		converg = zeros(Float64, maxState)
	end

	instance = zeros(Float64, maxState)

	iters = [0]
	vals = [0.0]

	indication = true

	modulo = 10
	totalTimes = 300000
	write(STDERR, "$totalTimes")


	for beans in 1:totalTimes
		 #not sure how to start the chain, doesn't really matter but still
		t = time()
		
		#get graph adj list for state that we're on currenlty
		arr = G[a,b]

		#get a random end state for us to jump to
		i = fRand()
		i = ceil( i/probs[f(a,b, maxState)] )   #check this line of code
		i = convert(Int64, i)
		#arr = G[a,b]
		c = arr[i]
	
		#update the number of times we've been at a state
		instance[c] += 1
		converg[c] = 1.0 / instance[c]

		#update the convergence distribution each time we get to a state
		indicator = true
		for j in 1:maxState
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
			write(STDOUT, "All states reached at iteration $beans and at state groupings $a $b $c\n")
			indication = false
		end

		#update the new difference and add it to the x and y axis information
		d = norm(converg-station)
		t = time() - t
		if beans % modulo == 0
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


	for i in 1:maxState
		print("My: $(converg[i]) \tStationary: $(station[i])\n")
		if converg[i] == 1
			numOnes += 1
		end
	end
	println("Num ONes: $numOnes")
	println("original norm: $d")

	for i in 1:maxState
		if instance[i] == 0
			converg[i] = 0
		end
	end

	d = norm(converg-station)
	println("difference in norm: $d")

	n = size(iters, 1)
	iters = iters[2:n]
	vals = vals[2:n]
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


function simuJump()

	inputFile = "$networkName-In.txt"
	probsFile = "$networkName-Probs.txt"
	G, M, maxState, probs = getMatrix(inputFile, probsFile)

	#getting information about the size of the eigenvectors

	n = size(M,1)

	station = readdlm("$networkName-StationaryDist.out")

	write(STDERR, "beginning simulation\n")
	@time iters, vals = SimulateSOMC(G, station, maxState, probs)

	plotInfo = [iters vals]
	writedlm("$networkName-Plot.txt", plotInfo)


	@time iters, vals = SimulateSOMC(G, station, maxState, probs, true)

	plotInfo = [iters vals]
	writedlm("$networkName-Plot1.txt", plotInfo)


	#plot function after all is done
	plot(iters, vals)
	show()
end

function printGraph()

	inputFile = "$networkName-In.txt"
	probsFile = "$networkName-Probs.txt"
	G, M, maxState, probs = getMatrix(inputFile, probsFile)
	println("$G")
	quit()
end

function main()

	inputFile = "$networkName-In.txt"
	probsFile = "$networkName-Probs.txt"
	G, M, maxState, probs = getMatrix(inputFile, probsFile)

	#getting information about the size of the eigenvectors

	n = size(M,1)

	colStoch = true#testColumnStochastic(M, n)

	if !colStoch
		println("Error, matrix is not column stochastic")
		quit()
	end
	#calculated, isolate, and print the eigenvalues and eigenvector information


	c, eigenSize = extractNormalizedEigenvector(M)
	writedlm("$networkName-EgVec.out", c)

	#calculating the one state stationary distribution, instead of over pairs of states
	station = extractStationaryDist(maxState, c, eigenSize)
	writedlm("$networkName-StationaryDist.out", station)

	#begin simulation, and then plot the convergence rate and print axis information to file

	write(STDERR, "beginning simulation\n")
	@time iters, vals = SimulateSOMC(G, station, maxState, probs)

	plotInfo = [iters vals]
	writedlm("$networkName-Plot0.txt", plotInfo)

	@time iters, vals = SimulateSOMC(G, station, maxState, probs, true)

	plotInfo = [iters vals]
	writedlm("$networkName-Plot1.txt", plotInfo)

	#plot function after all is done
	plot(iters, vals)
	show()
end

write(STDERR, "begin\n")
if networkName == "NAN"
	println("Change network name")
	quit()
end
main()
