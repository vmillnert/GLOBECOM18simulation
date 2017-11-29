#=

Author: Victor Millnert
email:  victor@control.lth.se

Main method for the simulations.


This one is tailored to investigate how lambda affect the utility, amount
of discarded packets, and amount of overallocation

=#
@everywhere include("./simulation.jl")
@everywhere include("./simulation_sota.jl")


@everywhere using LightGraphs
@everywhere using StatsBase
@everywhere using Distributions
@everywhere using Convex
@everywhere using SCS


include("./write_output.jl")


DAS_low = 0.6 # lower threshold  => scale down if below
DAS_high = 0.8 # upper threshold => scale up if above
DOA_ratio = 1.2 # input rate x this

# Let's vary lambda
lambda = 0.5


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

# The bounds for the randomly generated deadlines for each local nodes
# TODO: change this to instead represent the bounds for the deadlines for the paths
deadline_low  = 1.0  # time-units (in seconds)
deadline_high = 3.0  # time-units

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


M = 2

U_mean = SharedArray{Float64}(M)
A_mean = SharedArray{Float64}(M)
E_mean = SharedArray{Float64}(M)
Discarded = SharedArray{Float64}(M)
Overallocation = SharedArray{Float64}(M)

U_mean_DAS = SharedArray{Float64}(M)
A_mean_DAS = SharedArray{Float64}(M)
E_mean_DAS = SharedArray{Float64}(M)
Discarded_DAS = SharedArray{Float64}(M)
Overallocation_DAS = SharedArray{Float64}(M)

U_mean_DAS_AC = SharedArray{Float64}(M)
A_mean_DAS_AC = SharedArray{Float64}(M)
E_mean_DAS_AC = SharedArray{Float64}(M)
Discarded_DAS_AC = SharedArray{Float64}(M)
Overallocation_DAS_AC = SharedArray{Float64}(M)


U_mean_DOA = SharedArray{Float64}(M)
A_mean_DOA = SharedArray{Float64}(M)
E_mean_DOA = SharedArray{Float64}(M)
Discarded_DOA = SharedArray{Float64}(M)
Overallocation_DOA = SharedArray{Float64}(M)

U_mean_DOA_AC = SharedArray{Float64}(M)
A_mean_DOA_AC = SharedArray{Float64}(M)
E_mean_DOA_AC = SharedArray{Float64}(M)
Discarded_DOA_AC = SharedArray{Float64}(M)
Overallocation_DOA_AC = SharedArray{Float64}(M)


# let's do it as a short Monte Carlo simulation
@sync @parallel for m in 1:M

    # simulate our method
    U_mean[m], A_mean[m],
    E_mean[m], Discarded[m],
    Overallocation[m] = sim_graph(sim_time,
                                  P,SUNET, lambda,
                                  dt, Tend,t, N, peak,
                                  deadline_low, deadline_high,
                                  timeoverhead_low, timeoverhead_high,
                                  nominal_service_low, nominal_service_high,
                                  xi_low, xi_high);


    # simulate dynamic auto scaling without admission control
    DAS   = true
    DOA   = false
    AC    = false

    U_mean_DAS[m], A_mean_DAS[m],
    E_mean_DAS[m], Discarded_DAS[m],
    Overallocation_DAS[m] = sim_graph_sota(sim_time,
                                           P,
                                           DAS, DOA, AC,
                                           DAS_low, DAS_high,
                                           DOA_ratio,
                                           SUNET, lambda,
                                           dt, Tend,t, N, peak,
                                           deadline_low, deadline_high,
                                           timeoverhead_low, timeoverhead_high,
                                           nominal_service_low, nominal_service_high,
                                           xi_low, xi_high);

    # simulate dynamic auto scaling with admission control
    DAS   = true
    DOA   = false
    AC    = true

    U_mean_DAS_AC[m], A_mean_DAS_AC[m],
    E_mean_DAS_AC[m], Discarded_DAS_AC[m],
    Overallocation_DAS_AC[m] = sim_graph_sota(sim_time,
                                              P,
                                              DAS, DOA, AC,
                                              DAS_low, DAS_high,
                                              DOA_ratio,
                                              SUNET, lambda,
                                              dt, Tend,t, N, peak,
                                              deadline_low, deadline_high,
                                              timeoverhead_low, timeoverhead_high,
                                              nominal_service_low, nominal_service_high,
                                              xi_low, xi_high);

    # simulate dynamic auto scaling without admission control
    DAS   = false
    DOA   = true
    AC    = false

    U_mean_DOA[m], A_mean_DOA[m],
    E_mean_DOA[m], Discarded_DOA[m],
    Overallocation_DOA[m] = sim_graph_sota(sim_time,
                                           P,
                                           DAS, DOA, AC,
                                           DAS_low, DAS_high,
                                           DOA_ratio,
                                           SUNET, lambda,
                                           dt, Tend,t, N, peak,
                                           deadline_low, deadline_high,
                                           timeoverhead_low, timeoverhead_high,
                                           nominal_service_low, nominal_service_high,
                                           xi_low, xi_high);
    
    # simulate dynamic auto scaling with admission control
    DAS   = false
    DOA   = true
    AC    = true

    U_mean_DOA_AC[m], A_mean_DOA_AC[m],
    E_mean_DOA_AC[m], Discarded_DOA_AC[m],
    Overallocation_DOA_AC[m] = sim_graph_sota(sim_time,
                                              P,
                                              DAS, DOA, AC,
                                              DAS_low, DAS_high,
                                              DOA_ratio,
                                              SUNET, lambda,
                                              dt, Tend,t, N, peak,
                                              deadline_low, deadline_high,
                                              timeoverhead_low, timeoverhead_high,
                                              nominal_service_low, nominal_service_high,
                                              xi_low, xi_high);

    
    println("---------------------------------------")
    println("Iteration $m (of $M) done")
    println("---------------------------------------")

end # for m

u_mean = mean(U_mean)
a_mean = mean(A_mean)
e_mean = mean(E_mean)
discarded = mean(Discarded)
overallocation = mean(Overallocation)

u_var  = var(U_mean)
a_var  = var(A_mean)
e_var  = var(E_mean)
discarded_var  = var(Discarded)
overallocation_var  = var(Overallocation)


u_mean_DAS = mean(U_mean_DAS)
a_mean_DAS = mean(A_mean_DAS)
e_mean_DAS = mean(E_mean_DAS)
discarded_DAS = mean(Discarded_DAS)
overallocation_DAS = mean(Overallocation_DAS)

u_var_DAS  = var(U_mean_DAS)
a_var_DAS  = var(A_mean_DAS)
e_var_DAS  = var(E_mean_DAS)
discarded_var_DAS  = var(Discarded_DAS)
overallocation_var_DAS  = var(Overallocation_DAS)


u_mean_DAS_AC = mean(U_mean_DAS_AC)
a_mean_DAS_AC = mean(A_mean_DAS_AC)
e_mean_DAS_AC = mean(E_mean_DAS_AC)
discarded_DAS_AC = mean(Discarded_DAS_AC)
overallocation_DAS_AC = mean(Overallocation_DAS_AC)

u_var_DAS_AC  = var(U_mean_DAS_AC)
a_var_DAS_AC  = var(A_mean_DAS_AC)
e_var_DAS_AC  = var(E_mean_DAS_AC)
discarded_var_DAS_AC  = var(Discarded_DAS_AC)
overallocation_var_DAS_AC  = var(Overallocation_DAS_AC)


u_mean_DOA = mean(U_mean_DOA)
a_mean_DOA = mean(A_mean_DOA)
e_mean_DOA = mean(E_mean_DOA)
discarded_DOA = mean(Discarded_DOA)
overallocation_DOA = mean(Overallocation_DOA)

u_var_DOA  = var(U_mean_DOA)
a_var_DOA  = var(A_mean_DOA)
e_var_DOA  = var(E_mean_DOA)
discarded_var_DOA  = var(Discarded_DOA)
overallocation_var_DOA  = var(Overallocation_DOA)

u_mean_DOA_AC = mean(U_mean_DOA_AC)
a_mean_DOA_AC = mean(A_mean_DOA_AC)
e_mean_DOA_AC = mean(E_mean_DOA_AC)
discarded_DOA_AC = mean(Discarded_DOA_AC)
overallocation_DOA_AC = mean(Overallocation_DOA_AC)

u_var_DOA_AC  = var(U_mean_DOA_AC)
a_var_DOA_AC  = var(A_mean_DOA_AC)
e_var_DOA_AC  = var(E_mean_DOA_AC)
discarded_var_DOA_AC  = var(Discarded_DOA_AC)
overallocation_var_DOA_AC  = var(Overallocation_DOA_AC)


println("-----------------------")
println("      our method       ")
println("Utility mean: $(u_mean)")
println("Availability mean: $(a_mean))")
println("Efficiency mean: $(e_mean)")
println("% discarded: $(discarded)")
println("% overallocation: $(overallocation)")

println("-----------------------")
println("          DAS          ")
println("Utility mean: $(u_mean_DAS)")
println("Availability mean: $(a_mean_DAS))")
println("Efficiency mean: $(e_mean_DAS)")
println("% discarded: $(discarded_DAS)")
println("% overallocation: $(overallocation_DAS)")

println("-----------------------")
println("          DOA          ")
println("Utility mean: $(u_mean_DOA)")
println("Availability mean: $(a_mean_DOA))")
println("Efficiency mean: $(e_mean_DOA)")
println("% discarded: $(discarded_DOA)")
println("% overallocation: $(overallocation_DOA)")

println("-----------------------")
println("       DAS + AC        ")
println("Utility mean: $(u_mean_DAS_AC)")
println("Availability mean: $(a_mean_DAS_AC))")
println("Efficiency mean: $(e_mean_DAS_AC)")
println("% discarded: $(discarded_DAS_AC)")
println("% overallocation: $(overallocation_DAS_AC)")

println("-----------------------")
println("       DOA + AC        ")
println("Utility mean: $(u_mean_DOA_AC)")
println("Availability mean: $(a_mean_DOA_AC))")
println("Efficiency mean: $(e_mean_DOA_AC)")
println("% discarded: $(discarded_DOA_AC)")
println("% overallocation: $(overallocation_DOA_AC)")

save_comparison_sim_data(u_mean, e_mean, a_mean, discarded, overallocation,
                         u_var, e_var, a_var, discarded_var, overallocation_var,
                         u_mean_DAS, e_mean_DAS, a_mean_DAS, discarded_DAS, overallocation_DAS,
                         u_var_DAS, e_var_DAS, a_var_DAS, discarded_var_DAS, overallocation_var_DAS,
                         u_mean_DAS_AC, e_mean_DAS_AC, a_mean_DAS_AC, discarded_DAS_AC, overallocation_DAS_AC,
                         u_var_DAS_AC, e_var_DAS_AC, a_var_DAS_AC, discarded_var_DAS_AC, overallocation_var_DAS_AC,
                         u_mean_DOA, e_mean_DOA, a_mean_DOA, discarded_DOA, overallocation_DOA,
                         u_var_DOA, e_var_DOA, a_var_DOA, discarded_var_DOA, overallocation_var_DOA,
                         u_mean_DOA_AC, e_mean_DOA_AC, a_mean_DOA_AC, discarded_DOA_AC, overallocation_DOA_AC,
                         u_var_DOA_AC, e_var_DOA_AC, a_var_DOA_AC, discarded_var_DOA_AC, overallocation_var_DOA_AC)
