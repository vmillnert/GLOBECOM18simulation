#=

Author: Victor Millnert
email:  victor@control.lth.se


Just a simple script to try out the deadline splitting scheme using
optimization.

=#

using Convex, SCS

n = 3
eps = 1e-3

Pmax = [10; 35; 32]

# define the path-node matrix (element p,v is 1 if path p goes through
# node v)
PN = [2 2 1;
      1 0 1;
      0 1 1]

l = length(Pmax)

Delta = [100 90 90]'


x = Variable(n);
problem = maximize(sum(x./Delta)/n + minimum(x./Delta), eps < x, PN * x <= Pmax)
solve!(problem, SCSSolver(verbose=0))


println("--------------")
println("Problem 1: minimize average ratio + minimum ratio")
println("average slack: $(sum((Pmax .- PN*x.value)./Pmax)/l)")
println("average ratio: $(mean(x.value./Delta))")
println("minimum ratio: $(minimum(x.value./Delta))")
println("deadlines: $(x.value)")
println("--------------")



x2 = Variable(n);
problem2 = minimize(sum(Delta./x2)/n + maximum(Delta./x2), eps < x2, PN * x2 <= Pmax)
solve!(problem2, SCSSolver(verbose=0))

println("--------------")
println("Problem 2: minimize average slack")
println("average slack: $(sum((Pmax .- PN*x2.value)./Pmax)/l)")
println("average ratio: $(mean(x2.value./Delta))")
println("minimum ratio: $(minimum(x2.value./Delta))")
println("deadlines: $(x2.value)")
println("--------------")


# x2 = Variable(n);
# problem2 = minimize(sum(Pmax .- PN*x2), eps < x2, PN * x2 <= Pmax)
# solve!(problem2, SCSSolver(verbose=0))

# println("--------------")
# println("Problem 2: minimize average slack")
# println("average slack: $(sum((Pmax .- PN*x2.value)./Pmax)/l)")
# println("average ratio: $(mean(x2.value./Delta))")
# println("minimum ratio: $(minimum(x2.value./Delta))")
# println("deadlines: $(x2.value)")
# println("--------------")
