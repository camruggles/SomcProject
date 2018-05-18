


using PyPlot			#inclusion of relevant library

#the value of network name needs to be set accordingly for input and output files to be used properly

networkName = ARGS[1]		#the name of the network used for all the file naming ,should be NAN unless it's a copy fo the source
stateA = 621			#set to a dummy value
stateB = 1478			#set to a dummy value

totalTimes = parse(Int, ARGS[2])		#the number of iterations in the simulation

#=
	Sums over pairs of states to create a stationary distribution of normal size
	args: maxState: the largest state number or the number of states
		c: the eigenvector normalized
		eigenSize: the size of the eigenvector
	output: the stationary distribution
	@author Cameron RUggles
=#
function extractStationaryDist(maxState, c, eigenSize)
	#Sum over pairs of states
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

#=
	extracts the appropriate eigenvector that goes with the stationary distribution, e.g. the eigenvector with eigenvalue 1
	args: a sparse chain matrix
	output: the eigenvector and the size
	@author Cameron Ruggles

=#
function extractEigenvector(M)

	#calculated, isolate, and print the eigenvalues and eigenvector information
	println("Eigenvector calculation")
	@time eigenvalues,eigenvectors = eigs(M, nev=5)

	#find the eigenvector with eigenvalue of one
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

	#extract the correct eigenvector
	eigenSize = size(eigenvectors, 1)
	c = eigenvectors[1:eigenSize, eigen_id]

	println("Eigenvalues")
	println("$eigenvalues")

	println("Eigenvector Size")
	println("$eigenSize")

	c = convert(Array{Float64, 1}, c)
	return c, eigenSize

end

#=
	Extracts and eigenvector from a sparse matrix and then normalizes it so it represents the stationary distribution
	args: a sparse chain matrix
	output: the stationary distribution and the size of the eigenvector
	@author Cameron Ruggles

=#
function extractNormalizedEigenvector(M)

	#extract eigenvector as normal
	c, eigenSize = extractEigenvector(M)

	#find the sum of the eigenvector for normalization using matrix method
	#this will work better than manually summing if julia is parallelized for that
	ei = ones(1, eigenSize)
	inner = ei*c

	println("inner product : $inner")
	s = inner[1]
	c = c/s

	return c, eigenSize
end


#=
	This will take two input files, read in all the state information and convert states into state pairs for a SOMC
	args: input file for state combinations, and another file for transition probabilities, which are all the same in one column and sum to one.
	output: G : the graph adjacency list, needed for the simulation later
			M : the sparse chain matrix
			maxState : the maximum number of states in the chain, or the number of the largest state
			probs : a probability transition dictionary for the transition prob at each state pair
	@author Cameron Ruggles

=#
function getMatrix(inputFile, probsFile)

	#taking the input states
	A = readdlm(inputFile, Int64)
	maxState = maximum(A) + 1
	println("Max State: $maxState")

	#reading in probabilities
	probs = Dict{Int64, Float64}()

	#adding probabilities to the dictionary
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



	#exstracting a pair of states to be used later during the beginning of the simulation
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

#=
	This will test to see if a matrix is column stochastic
	args: matrix 'M', and the size of the matrix 'n'
	output: a boolean value true or false indiciating column stochasticity
	@author Cameron Ruggles

=#
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


#=
	This will simulate a Markov Chain's random walk (not sure if that's right)
	args: graph G of markov chain adjList, the stationary distribution, the maximum number of states or the number of the largest state,
			a dictionary of the transition probabilities from each state, typically equal
			and a boolean value, set to one to use a vector of ones during norm measurements
			leave as is to use a regular vector of zeros
	output: two column array, with x axis and y axis information needed to generate a plot discussing the convergence
		will also print to file the simulated distribution after n iterations,
		and other information about when each state has been simulated at least once
	@author Cameron Ruggles

	this is a markov chain simulation, the chain makes random jumps to different states and we 
	record the information to see how quickly it becomes accurate to the stationary distribution

=#
function SimulateSOMC(G, station, maxState, probs, oneVec = false)
	counter = 0.0  # this is used to determine how many transitions have been made
	∆ = 0 #this is never used, but can be used to stop the simulation once the target vector is similar enough to the stationary distribution
	a,b = stateA, stateB #arbitrary starting point

	#initialize a target vector to be used for measuring convergence
	converg = 0
	if oneVec
		converg = ones(Float64, maxState)
	else
		converg = zeros(Float64, maxState)
	end

	#initializing variables
	instance = zeros(Float64, maxState)

	iters = [0]
	vals = [0.0]

	indication = true

	modulo = 10

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
		converg[c] = instance[c] / counter

		#update the convergence distribution each time we get to a state
		indicator = true

		#=
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
		=#

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
		#uncomment for monitoring purposes
		#print("My: $(converg[i]) \tStationary: $(station[i])\n")
		if instance[i] == 0
			numOnes += 1
		end
	end
	println("Num ONes: $numOnes")

	#getting rid of the zeros used to initialize the vectors
	#makes the plots look better
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

#=
This is a copy of the main function
If the stationary distribution is known from a file
then we read it in and then just run the simulation
used for hammering out details with simulation
Eigenvector calculation is very expensive and avoided in this case

@author Cameron RUggles
args: no input other than networkname and it's associated files, as well as the stationary distribution file
output: writes to file the Plot information for zero and one simulation
=#

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
	writedlm("$networkName-Plot0.out", plotInfo)


	@time iters, vals = SimulateSOMC(G, station, maxState, probs, true)

	plotInfo = [iters vals]
	writedlm("$networkName-Plot1.out", plotInfo)


	#plot function after all is done
	plot(iters, vals)
	show()
end


#=
Main execution of the algorithm
args: input includes files used for matrix construction, based on the network name
output: writes to file the following information: eigenvectors, stationary distribution, zero and one convergence plots
@author Cameron Ruggles

creates matrix for SOMC
extracts eigenvector and stationary distribution
and then runs converging simulation
=#
function main()

	#matrix construction from files
	inputFile = "$networkName-In.txt"
	probsFile = "$networkName-Probs.txt"
	G, M, maxState, probs = getMatrix(inputFile, probsFile)

	#getting information about the size of the eigenvectors
	n = size(M,1)

	#testing to make sure the matrix is column stochastic
	colStoch = true #testColumnStochastic(M, n)

	if !colStoch
		println("Error, matrix is not column stochastic")
		quit()
	end


	#calculated, isolate, and print the eigenvalues and eigenvector information
	c, eigenSize = extractNormalizedEigenvector(M)

	#include if you want, this file is generally very large and once the eigenvector is extracted, the stationary distribution is a quick find
	#writedlm("$networkName-EgVec.out", c)

	#calculating the one state stationary distribution, instead of over pairs of states
	station = extractStationaryDist(maxState, c, eigenSize)
	writedlm("$networkName-StationaryDist.out", station)

	#begin simulation, and then plot the convergence rate and print axis information to file
	write(STDERR, "beginning simulation\n")

	#first simulation
	@time iters, vals = SimulateSOMC(G, station, maxState, probs)
	plotInfo = [iters vals]
	writedlm("$networkName-Plot0.out", plotInfo)


	#second simulation with a different metric for convergence
	@time iters, vals = SimulateSOMC(G, station, maxState, probs, true)
	plotInfo = [iters vals]
	writedlm("$networkName-Plot1.out", plotInfo)

	#plot function after all is done
	#plot(iters, vals)
	#show()
end


#=
The actual starting point of the code execution
Used to make sure printing to stderr is working
and to prevent myself from making the same stupid mistake of not renaming the global variable "networkName" over and over again
=#
write(STDERR, "begin\n")	#stderr check

if networkName == "NAN"		#stupidity prevention
	println("Change network name")
	quit()
end

main()						#function to be run, not always main
