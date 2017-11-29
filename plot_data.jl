
################################################
# ---------------------------------------------#
# ---------------------------------------------#
#          Plot the simulation results         #
# ---------------------------------------------# 
# ---------------------------------------------#
################################################



####
# Plot metrics as a function of lambda
close("all")

#=
# ----
# Plot the path-wise input rates
figure("Lambda vs. utility", figsize=(10,5))
PyPlot.plot(lam, mean(U_mean, 2), label="mean utility")
# PyPlot.plot(lam, Discarded, label="% discarded")
# PyPlot.plot(lam, Overallocation, label="% overallocation")
xlabel("lambda")
ylabel("Utility")
legend(loc="lower right")
title("Lambda vs. utility")
# ----    

# ----
# Plot the path-wise input rates
figure("Lambda vs. discarded packets", figsize=(10,5))
# PyPlot.plot(lam, U_mean, label="mean utility")
PyPlot.plot(lam, mean(Discarded, 2)*100, label="% discarded")
PyPlot.plot(lam, mean(Overallocation, 2)*100, label="% overallocation")
xlabel("lambda")
ylabel("percent")
legend(loc="lower right")
title("Lambda vs. discarded packets")
# ----    
=#

####
# Plot metrics as a function of time_overhead vs. deadline ratio
close("all")


# ----
# Plot the path-wise input rates
figure("Time overhead / deadline ration -- utility ", figsize=(10,5))
PyPlot.plot(dead_timeOH_ratio, mean(U_mean, 2), label="mean utility")
PyPlot.plot(dead_timeOH_ratio, mean(a_mean, 2), label="mean availability")
PyPlot.plot(dead_timeOH_ratio, mean(e_mean, 2), label="mean efficiency")
# PyPlot.plot(lam, Discarded, label="% discarded")
# PyPlot.plot(lam, Overallocation, label="% overallocation")
xlabel("ratio timeoverhead / deadline")
ylabel("Utility")
legend(loc="lower right")
title("Time overhead / deadline ration -- utility")
# ----    

# ----
# Plot the path-wise input rates
figure("Time overhead / deadline ratio", figsize=(10,5))
# PyPlot.plot(lam, U_mean, label="mean utility")
PyPlot.plot(dead_timeOH_ratio, mean(Discarded, 2)*100, label="% discarded")
PyPlot.plot(dead_timeOH_ratio, mean(Overallocation, 2)*100, label="% overallocation")
xlabel("ratio timeoverhead / deadline")
ylabel("percent")
legend(loc="lower right")
title("Time overhead / deadline ratio")
# ----    


#=



################################################
# ---------------------------------------------#
# ---------------------------------------------#
#          Plot the simulation results         #
# ---------------------------------------------# 
# ---------------------------------------------#
################################################
close("all")


# ----
# Plot the path-wise input rates
for p in 1:length(P)
    figure("Path-wise input rates", figsize=(10,5))
    PyPlot.plot(t, R[p,:]./dt, label="path $p: $(P[p])")
    xlabel("time (s)")
    ylabel("packets per second")
    legend(loc="lower right")
    title("Path-wise input rates")
end    
# ----    

# ----
# Plot the ref-signal
for v in 1:n
    figure("Ref-signal for number of instances", figsize=(10,5))
    PyPlot.plot(t, mref[v,:], label="mref_$v")
    xlabel("time (s)")
    ylabel("# instances")
    legend(loc="lower right")
    title("Ref-signal for number of instances")
end    
# ----    


# ----
# Plot the number of running machines
for v in 1:n
    figure("Number of instances", figsize=(10,5))
    PyPlot.plot(t, m[v,:], label="m_$v")
    xlabel("time (s)")
    ylabel("# instances")
    legend(loc="lower right")
    title("Number of instances")
end    
# ----    


# ----
# Plot the machine uncertainties
for v in 1:n
    figure("Average machine uncertainty (% of nominal service rate)", figsize=(10,5))
    PyPlot.plot(t, xi[v,:]./ns[v]*100, label="xi_$v(t)")
    xlabel("time (s)")
    ylabel("packets per second")
    legend(loc="lower right")
    title("Average machine uncertainty (% of nominal service rate)")
end    
# ----    


# # ----
# # Plot the routing weights over time
# for w in 1:l
#     figure("Routing weights for the edges", figsize=(10,5))
#     PyPlot.plot(t, W[w,:], label="edge $(E[w])")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Routing weights for the edges")
# end    
# # ----    

# # Plot the predicted routing weights over time
# for w in 1:l
#     figure("Routing weights for the edges", figsize=(10,5))
#     PyPlot.plot(t, What[w,:], label="what_$(E[w])")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Routing weights for the edges")
# end    
# # ----    

# # ----
# # Plot the external inputs
# for v in 1:n
#     figure("External inputs", figsize=(10,5))
#     PyPlot.plot(t, re[v,:], label="re_$v")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("External input")
# end    
# # ----    

# ----
# Plot the predicted external inputs
for v in 1:n
    figure("External inputs", figsize=(10,5))
    PyPlot.plot(t, re[v,:], label="re_$v")
    PyPlot.plot(t[Delta[v]+1:end], rehat[v,1:end-Delta[v]], label="rehat_$v")
    xlabel("time (s)")
    ylabel("packets per second")
    legend(loc="lower right")
    title("External input")
end    
# ----    


# # ----    
# # Plot the internal inputs
# for v in 1:n
#     figure("Internal input", figsize=(10,5))
#     PyPlot.plot(t, ri[v,:], label="ri_$v(t)")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Internal input")
# end    
# # ----

# ----    
# Plot the predicted internal inputs
for v in 1:n
    figure("Predicted internal input", figsize=(10,5))
    PyPlot.plot(t, ri[v,:], label="ri_$v(t)")
    PyPlot.plot(t[Delta[v]+1:end], rihat[v,1:end-Delta[v]], label="rihat_$v(t)")
    xlabel("time (s)")
    ylabel("packets per second")
    legend(loc="lower right")
    title("Predicted internal input")
end    
# ----    


# ----    
# Plot the service rate, arrival rate, and addmitance rate
for v in 1:n
    figure("Service, arrival, and admittance rate - F_$v", figsize=(10,5))
    PyPlot.plot(t, s[v,:]./dt, label="s_$v(t)")
    PyPlot.plot(t, a[v,:]./dt, label="a_$v(t)")
    PyPlot.plot(t, r[v,:]./dt, label="r_$v(t)")
    xlabel("time (s)")
    ylabel("packets per second")
    legend(loc="lower right")
    title("Service, arrival, and admittance rate - F_$v")
end    
# ----    

# # ----    
# # Plot the service rate vs the queue-size
# for v in 1:n
#     figure("Service rate vs queue size - F_$v", figsize=(10,5))
#     PyPlot.plot(t, scap[v,:]./dt, label="scap_$v(t)")
#     PyPlot.plot(t, A[v,:] - S[v,:], label="q_$v(t)")
#     PyPlot.plot(t, a[v,:]./dt, label="a_$v(t)")
#     PyPlot.plot(t, s[v,:]./dt, label="s_$v(t)")

#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Service rate vs queue size - F_$v")
# end    
# # ----    
 
# # ----    
# # Plot admittance and service curves
# for v in 1:n
#     figure("Admittance and service curves - F_$v", figsize=(10,5))
#     PyPlot.plot(t, S[v,:], label="S_$v(t)")
#     PyPlot.plot(t, A[v,:], label="A_$v(t)")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Admittance and service curves - F_$v")
# end    
# # ----    



## ## compute the admission policy difference
# Tmp = zeros(Float64, n, N)
# for v in 1:n
#     for i in 1:N-D[v]
#         Tmp[v,i] = A[v,i]+a[v,i] -( S[v,i] + Slb[v,i+D[v]] - Slb[v,i])
#     end
# end
    
# # ----    
# # Plot admittance and service curves
# for v in 1:n
#     figure("Admittance difference ", figsize=(10,5))
#     PyPlot.plot(t, Tmp[v,:], label="function $v")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Admittance difference")
# end    
# # ----  


# ### compute whether the packets meet their deadlines or not
# Tmp = zeros(Float64, n, N)
# for v in 1:n
#     for i in 1+D[v]:N
#         Tmp[v,i] = S[v,i] - A[v,i-D[v]]
#     end
# end
# # ----    
# # Plot deadline difference curve
# for v in 1:n
#     figure("Deadline difference curve (+ meet deadline, - miss) ", figsize=(10,5))
#     PyPlot.plot(t, Tmp[v,:], label="function $v")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Deadline difference curve (+ meet deadline, - miss)")
# end    
# # ----    


# # ----    
# # Plot utility functions
# for v in 1:n
#     figure("Utility function - F_$v ", figsize=(10,5))
#     PyPlot.plot(t, ua[v,:], label="ua_$v(t)")
#     PyPlot.plot(t, ue[v,:], label="ue_$v(t)")
#     PyPlot.plot(t, u[v,:], label="u_$v(t)")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Utility function - F_$v")
# end    
# # ----    

# ----    
# Plot average utility
figure("Average efficiency ", figsize=(10,5))
PyPlot.plot(t, Ue, label="Ue(t)")
xlabel("time (s)")
ylabel("packets per second")
legend(loc="lower right")
title("Average efficiency")

figure("Average availability ", figsize=(10,5))
PyPlot.plot(t, Ua, label="Ua(t)")
xlabel("time (s)")
ylabel("packets per second")
legend(loc="lower right")
title("Average availability")

figure("Average utility function ", figsize=(10,5))
PyPlot.plot(t, U, label="U(t)")
xlabel("time (s)")
ylabel("packets per second")
legend(loc="lower right")
title("Average utility function")
# ----    

=#
