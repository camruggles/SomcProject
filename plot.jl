

A = readdlm("correctConvRate.out")
n = size(A, 1)

x = A[1:n, 1]
y = A[1:n, 2]

using PyPlot
plot(x,y)
show()
