#=

Author: Victor Millnert
email:  victor@control.lth.se

This one is tailored to investigate how the ratio between the
time-overhead and the local node-deadlines affect the utility, amount
of discarded packets, and amount of overallocation

=#

@everywhere include("./simulation_deadline_ratio.jl")
include("./write_output.jl")

@everywhere using LightGraphs
@everywhere using StatsBase
@everywhere using Distributions


# Specify the paths
P = [(1,2,3,2,1);
     (2,3);
     (1,3)]


# Decide if we should use the SUNET data or not
SUNET = true


# Set up the length of the simulation-time - this is the data gathered
# from SUNET, i.e., the 1h data, 2h data, or the 6h data.
sim_time = 1
# sim_time = 2
# sim_time = 6



# Set up the time-parameters

dt = 1e-2 # simulation step-time
tmp_input, tmp_t = GetData.getInputData(sim_time, dt, 1e6, false, false)
Tend = round(tmp_t[end], round(Int64, 1/dt))
t = 0:dt:Tend
N = length(t)


# Peak of the input rate for the different packet flows. The
# aggregatet flow from multiple external input might of course be
# larger than this.
peak = 10.0e4*dt


# The bounds for the randomly generated time-overheads for each local
# node
timeoverhead_low  = 30.0 # time-units (in seconds)
timeoverhead_high = 2*60.0 # time-units


# The bounds for the randomly generated nominal service rate for each
# local node
nominal_service_low  = peak/5 # so we need more than 5 instances per input flow
nominal_service_high = peak/10 # so we need less than 10 instances per input flow

# The bounds for the randomly generated machine uncertainty bounds for each local
# node
xi_low  = -0.3 # minus 30% of the nominal service rate
xi_high = 0.3  # plus 30% of the nominal service rate


K = 2
M = 2

U_mean = SharedArray{Float64}(K, M)
A_mean = SharedArray{Float64}(K, M)
E_mean = SharedArray{Float64}(K, M)
Discarded = SharedArray{Float64}(K, M)
Overallocation = SharedArray{Float64}(K, M)

# Let's vary lambda
lambda = 0.5*ones(K)


# ####################
# Let's vary the ratio between local deadline and time-overhead
# deadline_ratio = time_overhead / local node deadlie
deadline_ratio = logspace(0,3,K)
# ###################



@sync @parallel for m in 1:M

    # let's do it as a short Monte Carlo simulation
    for k in 1:K

        
        # Let's vary lambda and see what happens!
        
        lam = lambda[k]        

        U_mean[k,m], A_mean[k,m],
        E_mean[k,m], Discarded[k,m],
        Overallocation[k,m] = sim_graph_deadline_ratio(sim_time,
                                        P,SUNET, lam,
                                        dt, Tend,t, N, peak,
                                        deadline_ratio[k],
                                        timeoverhead_low, timeoverhead_high,
                                        nominal_service_low, nominal_service_high,
                                        xi_low, xi_high);


        println("---------------------------------------")
        println("Iteration $k (of $K), $m (of $M) done")
        println("---------------------------------------")

    end
end


u_mean = mean(U_mean,2)
a_mean = mean(A_mean,2)
e_mean = mean(E_mean,2)
overallocation = mean(Overallocation,2)
discarded = mean(Discarded,2)

u_var  = var(U_mean,2)
a_var  = var(A_mean,2)
e_var  = var(E_mean,2)
discarded_var  = var(Discarded,2)
overallocation_var  = var(Overallocation,2)


println("Utility mean: $(u_mean)")
println("Availability mean: $(a_mean))")
println("Efficiency mean: $(e_mean)")
println("% discarded: $(discarded*100) %")
println("% overallocation: $(overallocation*100) %")


save_deadline_ratio_sim_data(deadline_ratio,
                             u_mean, e_mean, a_mean, discarded, overallocation,
                             u_var, e_var, a_var, discarded_var, overallocation_var)
