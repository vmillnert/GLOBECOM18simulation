#=

Author: Victor Millnert
email:  victor@control.lth.se


Simulation for evaluating the automatic service- and admission
controller developed for a NFV graph in the paper "Controlling NFV
Graphs with Hard Constraints".

This simulation is tailored towards a single long simulation and note
a Monte Carlo simulation. The difference is that this stores a lot of
data that can be used for plotting for instance the machine
uncertainties, service capacities, etc.

=#

include("./GetData.jl")
include("./write_output.jl")

using LightGraphs
using StatsBase
using Distributions
using PyPlot
using Convex
using SCS


# -------------------------------------------------
# Specify the user input
# -------------------------------------------------

P = [(1,2,3,2,1);
     (1,3);
     (2,3)]


# Decide if we should use the SUNET data or not
SUNET = true


# Set up the length of the simulation-time - this is the data gathered
# from SUNET, i.e., the 1h data, 2h data, or the 6h data.
# sim_time = 1
# sim_time = 2
sim_time = 6


###############
# Set up the time-parameters

dt = 1e-2 # simulation step-time
tmp_input, tmp_t = GetData.getInputData(sim_time, dt, 1e6, false, false)
Tend = round(tmp_t[end], round(Int64, 1/dt))
t = 0:dt:Tend
N = length(t)
###############


# Peak of the input rate for the different packet flows. The
# aggregatet flow from multiple external input might of course be
# larger than this.
peak = 10.0e4*dt

# The bounds for the randomly generated deadlines for each local nodes
# TODO: change this to instead represent the bounds for the deadlines for the paths
deadline_low  = 1.0  # time-units (in seconds)
deadline_high = 3.0  # time-units

# The bounds for the randomly generated time-overheads for each local
# node
time_overhead_low  = 30 # time-units (in seconds)
time_overhead_high = 2*60 # time-units


# The bounds for the randomly generated nominal service rate for each
# local node
nominal_service_low  = peak/5 # so we need more than 5 instances per input flow
nominal_service_high = peak/10 # so we need less than 10 instances per input flow

# The bounds for the randomly generated machine uncertainty bounds for each local
# node
xi_low  = -0.3 # minus 30% of the nominal service rate
xi_high = 0.3  # plus 30% of the nominal service rate


# Specify lambda for the functions
lambda = 0.5



# -------------------------------------------------
# Get the path-wise input rates from the SUNET data
# -------------------------------------------------

Rpath = zeros(Float64, length(P), N)
for p in 1:length(P)

    if SUNET
        ############################################################
        # If we want the SUNET data as input
        Rpath[p,:] = GetData.getInputData(sim_time, dt, peak, false, false)[1];
        ############################################################
    else
        
        ###############
        # If we want sinus-waves as input
        time_constant = rand(linspace(100,1000,100))
        time_shift = rand(linspace(-pi,pi,100))

        for i in 1:N
            # Get the flow weight from the SUNET data
            
            # without Gaussian noise
            Rpath[p,i] = peak*1/2 + cos(t[i]/time_constant + time_shift)*peak*1/3

            # This part is with some Gaussian noise
            # Rpath[p,i] = peak*1/2 + max(0,rand(Normal(abs(cos(t[i]/time_constant+
            #                                              time_shift))*peak*1/3,
            #                                      peak/1e2)))
        end
        ###################
    end
end



# Tend = 1e3 # stop time for the simulation
# t = 0:dt:Tend # simulation length (in seconds)
# N = length(t) # number of simulation steps


####################################
# Define the various attributes of the VNFs
dk = round.(Int64, 10/dt) # number of steps the derivative is
                                 # computed over



# ---------------------------------
# Set up the graph
# ---------------------------------

# Find the number of unique nodes
n = 0
for p in P
    for i in p
        n = max(n,i)
    end
end

# create the graph
g = DiGraph(n)

# add the edges
for p in P
    for i in 1:length(p)-1
        add_edge!(g, p[i],p[i+1])
    end
end

# number of links in the graph
l = ne(g)

# # get the adjacency matrix
# Adj = adjacency_matrix(g)


# Set up a node-link incidence matrix
# This one should be given by the iterable given from edges(g)
# NL is of size #nodes x #links
NL = spzeros(n, l)
E = Vector{Tuple{Int64,Int64}}(l) # a
# list
# containing
# all
# the
# edges
i = 1
for e in edges(g)
    E[i] = src(e),dst(e)
    NL[src(e),i] = 1
    NL[dst(e),i] = -1
    i = i+1
end

########
# Set up a link-path incidence matrix
# LP is of size #links x #paths

# It is important to note here that if a path goes over a link
# multiple times, that link should have a higher number!

LP = spzeros(l, length(P))
j = 1
for p in P # iterate through each path
    for i in 1:length(p)-1 # iterate through each link on the path

        # add each link on the path to the LP-matrix
        k::Int64 = find(x -> x==(p[i],p[i+1]), E)[1]            
        LP[k,j] = LP[k,j] + 1 # increment the k,j'th entry of the
        # link-path matrix by one
    end
    j = j+1
end


# ----------
# Compute the node-path incidence matrix

# PN is of size #paths x #paths. The PN-matrix is a place holder of
# how many times a path passes through a specific node

# We define one node-path matrix for the input flows (NPin) and one
# node-path matrix for the out flows (NPout)


PN = spzeros(Int, length(P), n)
NPin = spzeros(Int, n, length(P))
NPout = spzeros(Int, n, length(P))

j = 1
for p in P
    input_i::Int64 = p[1]
    output_i::Int64 = p[end]
    NPin[input_i, j] =  1
    NPout[output_i, j] =  -1
    for i in 1:length(p)
        PN[j,p[i]] = PN[j, p[i]] + 1
    end
    j = j+1
end

# @show full(PN)

# ---------





# ---------------------------------------------
# Set up the matrices needed for the simulation
# ---------------------------------------------


# The routing weights for the links in the graph.
#W[e,i] correspond to the routing weight for edge e at time-index i
W = zeros(Float64, l, N)
What = zeros(Float64, l, N)


# arrival rate matrices
re = zeros(Float64, n, N)
ri = zeros(Float64, n, N)
r  = zeros(Float64, n, N)
R  = zeros(Float64, n, N)

rehat = zeros(Float64, n, N)
rihat = zeros(Float64, n, N)
rhat  = zeros(Float64, n, N)



# admittance rate matrices
a = zeros(Float64, n, N)
A = zeros(Float64, n, N)

# service rate matrices
s    = zeros(Float64, n, N)
S    = zeros(Float64, n, N)
scap = zeros(Float64, n, N)
Slb  = zeros(Float64, n, N)
Scap = zeros(Float64, n, N)

ns = rand(linspace(nominal_service_low,
                                     nominal_service_high, 1e4), n) # draw
# a
# random
# nominal
# service
# rate

# machine matrices
m        = zeros(Int64, n, N)
mref     = zeros(Int64, n, N)
kappa    = zeros(Float64, n, N)

# machine uncertainties
xilb = rand(linspace(xi_low, 0, 1e4), n).*ns # draw a lower bound by
# random withing 0-30 %
# of the nominal service
# rate
xiub = rand(linspace(0, xi_high, 1e4), n).*ns # draw an upper bound at
# random from 0-30% of
# the nominal service
# rate
xi   = zeros(Float64, n, N)
xiseed = zeros(Float64, 100, n)

for v in 1:n
    xiseed[:,v] = linspace(xilb[v], xiub[v], 100)
end

# utility functions
ue  = zeros(Float64, n, N) # efficiency metric
ua  = zeros(Float64, n, N) # availability metric
u   = zeros(Float64, n, N) # utility function
lam = lambda*ones(Float64, n) # relative importance of
# efficiency and
# availability

Ue  = zeros(Float64, N) # average efficiency metric
Ua  = zeros(Float64, N) # average availability metric
U   = zeros(Float64, N) # average utility


# randomly generate the time-overheads
Delta_float = rand(linspace(time_overhead_low,
                                              time_overhead_high,
                                              1e4), n)


################################################
# ---------------------------------------------
# Solve the deadline splitting problem

# Generate random path-wise end-to-end deadlines
# one per path, randomly generated within the interval specified by
# [deadline_low, deadline_high]
Pmax = rand(linspace(deadline_low,
                     deadline_high,
                     1e4), length(P))

# specify the variables for the optimization problem
D_x = Variable(n);


D_eps = dt*3

# solve the problem
problem = maximize(sum(D_x./Delta_float)/n +
                   minimum(D_x./Delta_float),
                   D_eps < D_x, PN * D_x <= Pmax)

solve!(problem, SCSSolver(verbose=0))


# convert the Delta-matrix to be the number of time-steps instead of
# the actual time. This speeds up the simulation significantly
Delta = round.(Int64, Delta_float ./ dt)

# extract the node-deadlines from the solution of the cone-program
D = floor.(Int64, D_x.value[:,1] ./ dt)

# @show D
# @show D_x.value

# ---------------------------------------------
################################################


################################################
# ---------------------------------------------#
# ---------------------------------------------#
#              Run the simulation              #
# ---------------------------------------------# 
# ---------------------------------------------#
################################################

xi_tmp = 0

# set the initial ref-value to be given by the total input of the
# paths flowing through that node

temp = zeros(Float64, n)
for p in 1:length(P)
    for v in P[p]
        temp[v] = temp[v] + Rpath[p,1]
    end
end
for v in 1:n
    for i in 1:Delta[v]
        m[v,i] = round.(Int64, temp[v] / (ns[v] + xilb[v]))
        mref[v,i] = round.(Int64, temp[v] / (ns[v] + xilb[v]))
    end
    # set up some initial machine uncertainties
    xi_tmp = rand(Normal()).*(xiub[v]-xilb[v])
    xi[v,1] = min(max(xilb[v], xi_tmp), xiub[v])

end

# We have to compute Silb for the first time-units
for v in 1:n
    for i in 2:Delta[v]+1
        Slb[v,i] = Slb[v,i-1] + m[v,i]*(ns[v] + xilb[v])
    end
end


# ----------------------------------------------
# Compute the routing weights for the simulation
# ----------------------------------------------

println("Computing routing weights")
for i in 1:N
    
    if mod(i, round.(Int,N/4)) == 0
        println("$(round.(Int, i/N*100)) % done")
    end

    for v in 1:n
        # we're at node v. Iteration index i

        # we have to compute the tmp-value from all the outflowing edges as
        # well as the special external output
        tot_input::Float64 = 0.0

        # we start by going through all the links going out of node i
        # we check if the link is positive
        # we check what path it belongs to
        for e in find(x-> x==1, NL[v,:]) # this extract the outgoing links
            for p in find(x-> x>0, LP[e,:]) # this extracts the paths
                # this link belongs to

                # note here that an input flow might flow over the
                # same link multiple times
                tot_input = tot_input + Rpath[p,i]*LP[e,p]
            end
        end

        # we also have to consider the external outputs of node i
        for p in find(x-> x==-1, NPout[v,:]) # extracts the paths that have
            # external outputs from node i
            tot_input = tot_input + Rpath[p,i]
        end

        # we can now normalize the weights for the outgoing links
        for e in find(x-> x==1, NL[v,:]) # this extract the outgoing links
            for p in find(x-> x>0, LP[e,:]) # this extracts the path
                # this link belongs to

                # again, one should note that a path might flow over a
                # link multiple times
                W[e,i] = W[e,i] + LP[e,p]*Rpath[p,i]/tot_input
            end
        end
    end
end

# first predicted routing weight
What[:,1:4] = W[:,1:4]

# intialize the predicted internal and external inputs
init = 4
for v in 1:n

    # predict all internal input rate
    for e in find(x-> x==1, NL[v,:]) # find all incoming links    
        for p in find(x-> x>0, LP[e,:]) # find all paths flowing over this link

            rihat[v,1:init] = rihat[v,1:init] + LP[e,p]*Rpath[p,1:init]

        end
    end

    # predict all external input rates
    for p in find(x-> x==1, NPin[v,:]) # find all external inputs to this node
        rehat[v,1:init] = rehat[v,1:init] + Rpath[p,1:init]
    end
    rhat[v,1:init] = rehat[v,1:init] + rihat[v,1:init]

end

# Some temporary variables used later on
x = 0.0 
mdiff = 0 
imin = 1
ein = 0
eout = 0

println("Running simulation")
for i in 1:N

    if mod(i, round.(Int,N/4)) == 0
        println("$(round.(Int, i/N*100)) % done")
    end

    
    # iterate through each of the VNF nodes
    for v in 1:n

        ##
        #1. Gather the external inputs
        ##
        
        # this is easily found by the +1 in the node-path incidence matrix
        for p in find(x-> x==1, NPin[v,:]) # extracts the paths that
            # have external outputs
            # from node i
            re[v,i] = re[v,i] + Rpath[p,i]
        end

        ##
        #2. Compute the admittance rate (i.e., the admission control)
        ##
        r[v,i] = re[v,i] + ri[v,i] # compute total arrival rate
        a[v,i] = r[v,i] # we start by admitting packets
        if i + D[v] < N
            if A[v,i] + a[v,i] > S[v,i] + Slb[v,i+D[v]] - Slb[v,i]
                # there is a chance to miss deadlines => discard
                # incoming packets
                # a[v,i] = 0.0
                if i +D[v] - Delta[v] > 0
                    a[v,i] = min(r[v,i], mref[v,i+D[v]-Delta[v]]*(ns[v] + xilb[v]))
                else
                    a[v,i] = 0.0
                end
                
            end
        end

        ##
        #3. Compute the service rate
        ##

        # update the number of running instances
        if i > Delta[v]
            m[v,i] = mref[v,i-Delta[v]]
        end
        
        # First we need to see whether we should change the machine
        # uncertainty
        
        # if the number of running machines differ from last
        # time-step, we should update and simulate new machine
        # uncertainties

        if i > 1
            mdiff = abs(m[v,i] - m[v,i-1])

            # get some new random machine uncertainty
            # x = rand(Normal()).*(xiub[v]-xilb[v])

            x = rand(xiseed[:,v])

            # make the difference proportional to the change in the number of instances
            x = (xi[v,i-1]*m[v,i-1] + x*mdiff) / (mdiff + m[v,i-1]) 

            xi[v,i] = min(max(xilb[v], x), xiub[v])
        end

        x = rand(Normal(xi[v,i], 1/1000))
        xi[v,i] = min(max(xilb[v], x), xiub[v])
        
        # compute maximum service capacity, scap
        scap[v,i] = m[v,i]*(ns[v] + xi[v,i])

        
        # process the packets
        s[v,i] = min(scap[v,i], A[v,i] - S[v,i] + a[v,i])

        ##
        #4. Compute the ref-signal for the service rate
        ##        
        # compute the ref-signal
        if i > 1
            kappa[v,i] = rhat[v,i] / (ns[v] + xi[v,i])

            if (lam[v]*floor(kappa[v,i])*ceil(kappa[v,i]) +
                (1-2*lam[v])*kappa[v,i]*ceil(kappa[v,i]) >= (1-lam[v])*kappa[v,i]^2)

                mref[v,i] = floor(Int64, kappa[v,i])
            else
                # @show i, v, xi[v,i], kappa[v,i]
                mref[v,i] = ceil(Int64, kappa[v,i])
            end
        end        
        # always keep at least a single instance on
        mref[v,i] = max(mref[v,i], 1)


        ##
        #5. PREDICT future input rate, internal input rates, external
        #   input rates, and future routing weights
        ##

        # This should be done for the links going out from this node
        if i > 3 && i < N
            imin = max(1, i-dk)
            # predict future routing weights
            for e in find(x-> x>0, NL[v,:]) # get the outgoing links of this node
                
                What[e,i+1] = (W[e,i] + Delta[v] * (W[e,i] -
                                                    W[e,imin]) / (i + 1 - imin))

            end

            
            # predict the internal input rates
            for e in find(x-> x<0, NL[v,:]) # get all the incoming
                # links to this node

                ein = E[e][1] # the incoming neighbor of link e

                rihat[v,i+1] = (rihat[v,i+1] + What[e,i]*
                                min(mref[ein,i-1]*(ns[ein]+xi[ein,i-1]),
                                    rhat[ein,i-1]))
                
            end

            # predict the external input rates

            # ####
            # Remember! This node has no idea of what the different
            # incoming paths are. Therefore it is assumed that it
            # cannot differentiate between the different flows, but
            # only predict the aggregated external input rate.
            # ####
            
            rehat[v,i+1] = (re[v,i] + Delta[v]*(re[v,i] - re[v,imin]) / (i + 1 - imin))
            rhat[v,i+1] = rihat[v,i+1] + rehat[v,i+1]

        end
        


        if i<N
            ##
            #5. Compute the integrals
            ##
            A[v,i+1] = A[v,i] + a[v,i]
            S[v,i+1] = S[v,i] + s[v,i]
            R[v,i+1] = R[v,i] + r[v,i]
            Scap[v,i+1] = Scap[v,i] + scap[v,i]
            if i + Delta[v] < N
                Slb[v,i+1+Delta[v]] = Slb[v,i+Delta[v]] + mref[v,i]*(ns[v] + xilb[v])
            end


            ##
            #6. Send the internal outputs
            ##
            for e in find(x-> x>0, NL[v,:]) # get all the outgoing links from this node
                eout = E[e][2] # the outgoing neighbor of link e

                ri[eout,i+1] = ri[eout,i+1] + W[e,i]*s[v,i]
                
            end
            
            
        end

        ##
        #7. Compute the utility functions
        ##

        if r[v,i] > 0.0
            ua[v,i] = s[v,i]/r[v,i]
        end

# if we miss deadlines => 0 availability
if i > D[v] && S[v,i] < A[v,i-D[v]]
    ua[v,i] = 0.0
end

ue[v,i] = s[v,i]/scap[v,i]
# u[v,i]  = lam[v]*ua[v,i] + (1-lam[v])*ue[v,i]
u[v,i]  = ua[v,i]*ue[v,i]

end # for v in 1:n
U[i]  = sum(u[:,i])/n
Ue[i] = sum(ue[:,i])/n
Ua[i] = sum(ua[:,i])/n


end # for i in 1:N


# Compute % discarded packets
local disc = mean(1-A[:,end]./R[:,end])

# compute % overallocation
local overall = mean(Scap[:,end]./S[:,end]-1)

U_mean = mean(U,2)

Delta_max = maximum(Delta)

save_long_sim_data(Delta_max, t, Rpath, r, scap, m, xi, U_mean, a)
