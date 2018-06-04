


using PyPlot			#inclusion of relevant library

#the value of network name needs to be set accordingly for input and output files to be used properly

networkName = ARGS[1]		#the name of the network used for all the file naming
start_state_one = 621			#set to a dummy value
start_state_two = 1478			#set to a dummy value

totalIterations = parse(Int, ARGS[2])		#the number of iterations in the simulation

#=
	arguments:
		maxState : the index of the largest state
		c : the eigenvector
		eigenSize: the size of the eigenvector
	returns:
		the stationary distribution
=#
function extractStationaryDist(maxState, c, eigenSize)
	#Sum over pairs of states
	station = zeros(Float64, maxState, 1)
	for i in 1:maxState
		for j in i:maxState:eigenSize
			station[i] += c[j]
		end
	end

	return station
end

#=
	arguments:
		maxState : the index of the largest state
		c : the eigenvector
		eigenSize: the size of the eigenvector
	returns:
		the stationary distribution
=#
function extractRowStationaryDist(maxState, c, eigenSize)
	station = zeros(Float64, maxState, 1)
	for i in 1:eigenSize
		j = div(i-1, maxState) + 1
		station[j] += c[i]
	end

	return station
end

#=
	extracts the appropriate eigenvector that goes with the stationary distribution, e.g. the eigenvector with eigenvalue 1
	args: a sparse markov chain transition matrix
	output: the eigenvector and the size

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
=#
function extractNormalizedEigenvector(M)

	#extract eigenvector as normal
	c, eigenSize = extractEigenvector(M)

	#find the sum of the eigenvector for normalization using matrix method
	#this will work better than manually summing if julia is parallelized for that
	ei = ones(1, eigenSize)
	inner = ei*c

	println("inner product : $inner") #Should be one.
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
		i1 = statePairIndex(a,b,maxState)
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

	#exstracting a pair of states to be used later during the beginning of the simulation
	psize = size(A, 1)
	inl = psize * rand()
	inl = round(Int64, inl)

	global start_state_one = convert(Int64, A[inl,1]) + 1
	global start_state_two = convert(Int64, A[inl,2]) + 1

	#intaking input array and then calculating numbers for state pairs
	for i in 1:n
		#intaking the 3 variables
		a = convert(Int64, A[i,1])+1
		b = convert(Int64, A[i,2])+1
		c = convert(Int64, A[i,3])+1

		#converting the indices to states
		i1 = statePairIndex(a,b, maxState)
		i2 = statePairIndex(b,c, maxState)

		#recording the information in arrays for quicker Spare matrix construction
		arg1[i] = i1
		arg2[i] = i2
		arg3[i] = probs[i1]

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
	This will simulate a Markov Chain's random walk
	args: G: graph of markov chain adjList, 
			station: the stationary distribution, 
			maxState: the maximum number of states or the number of the largest state,
			probs: a dictionary of the transition probabilities from each state, typically equal
			oneVec: and a boolean value, set to one to use a vector of ones during norm measurements
			leave as is to use a regular vector of zeros
	output: two column array, with x axis and y axis information needed to generate a plot discussing the convergence
		will also print to file the simulated distribution after n iterations,
		and other information about when each state has been simulated at least once
=#
function SimulateSOMC(G, station, maxState, probs, oneVec = false)
    println("States: $start_state_one, $start_state_two")
	counter = 0.0  # this is used to determine how many transitions have been made
	a,b = start_state_one, start_state_two #arbitrary starting point

	#initialize a target vector to be used for measuring convergence
	simulatedDistribution = zeros(Float64, maxState)
	if oneVec
		simulatedDistribution = ones(Float64, maxState)
	end

	#initializing variables
	state_instances = zeros(Float64, maxState)
	unspoken = true #for printing information about a one time occurance

	iters = [0]
	vals = [0.0]

	printFrequency = 10 #use to determine how often to print some information


	for iterations in 1:totalIterations
		#get graph adj list for state that we're on currenlty
		arr = G[a,b]

		#get a random end state for us to jump to
		i = fRand()
		i = ceil( i/probs[statePairIndex(a,b, maxState)] )   #check this line of code
		i = convert(Int64, i)

		c = arr[i]
	
		#update the number of times we've been at a state
		state_instances[c] += 1
		simulatedDistribution[c] = state_instances[c] / counter

		

		#update the new difference and add it to the x and y axis information
		d = norm(simulatedDistribution-station)

		#print out information about the difference ever (modulo) transitions
		if iterations % printFrequency == 0
			append!(iters, iterations)
			append!(vals, d)
		end

		#print when the simulation falls below a certain threshold for vector difference
		#to deactivate, just set unspoken to false at the top of this function
		if d < 0.02 && unspoken
			println("Threshold: $iterations")
			unspoken = false
		end

		#update states for grabbing another random state during next iteration
		a=b
		b=c
		counter += 1
	end


	#clip the zeroes off the front of the two arrays to make the plots look better
	n = size(iters, 1)
	iters = iters[2:n]
	vals = vals[2:n]
	return iters, vals


end


#calculates the index of the state pairs
function statePairIndex(x, y, v)
	(x-1) * v + y
end


#Generates a number between 0 and 1, because I didn't know what rand() was
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
	writedlm("$networkName-Plot.out", plotInfo)

	#plot function after all is done
	plot(iters, vals)
	show()
end


#=
Main execution of the algorithm
args: input includes files used for matrix construction, based on the network name
output: writes to file the following information: eigenvectors, stationary distribution, zero and one convergence plots

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
	M = 0; gc()

	#calculating the one state stationary distribution, instead of over pairs of states
	station = extractRowStationaryDist(maxState, c, eigenSize)
	writedlm("$networkName-StationaryDist.out", station)

	#begin simulation, and then plot the convergence rate and print axis information to file
	write(STDERR, "\n\nbeginning simulation\n")

	#first simulation
	@time iters, vals = SimulateSOMC(G, station, maxState, probs)
	plotInfo = [iters vals]
	writedlm("$networkName-Plot.out", plotInfo)

	#plot function after all is done
	plot(iters, vals)
	show()
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
