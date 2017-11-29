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
    PyPlot.plot(t, Rpath[p,:]./dt, label="path $p: $(P[p])")
    xlabel("time (s)")
    ylabel("packets per second")
    legend(loc="lower right")
    title("Path-wise input rates")
end    
# ----    

# # ----
# # Plot the ref-signal
# for v in 1:n
#     figure("Ref-signal for number of instances", figsize=(10,5))
#     PyPlot.plot(t, mref[v,:], label="mref_$v")
#     xlabel("time (s)")
#     ylabel("# instances")
#     legend(loc="lower right")
#     title("Ref-signal for number of instances")
# end    
# # ----    


# ----
# Plot the input rates vs. service capacities
for v in 1:n
    figure("Input rate vs. capacity", figsize=(10,5))
    PyPlot.plot(t, r[v,:], label="r_$v(t)")
    PyPlot.plot(t, scap[v,:], label="smax_$v(t)")
    xlabel("time (s)")
    ylabel("packets per second")
    legend(loc="lower right")
    title("Input rate vs. capacity")
end    
# ----    

# ----
# Plot the number of running instances
for v in 1:n
    figure("Number of instances", figsize=(10,5))
    PyPlot.plot(t, m[v,:], label="m_$v")
    xlabel("time (s)")
    ylabel("# instances")
    legend(loc="lower right")
    title("Number of instances")
end    
# ----    


# # ----
# # Plot the machine uncertainties
# for v in 1:n
#     figure("Average machine uncertainty (% of nominal service rate)", figsize=(10,5))
#     PyPlot.plot(t, xi[v,:]./ns[v]*100, label="xi_$v(t)")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Average machine uncertainty (% of nominal service rate)")
# end    
# # ----    


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

# # ----
# # Plot the predicted external inputs
# for v in 1:n
#     figure("External inputs", figsize=(10,5))
#     PyPlot.plot(t, re[v,:], label="re_$v")
#     PyPlot.plot(t[Delta[v]+1:end], rehat[v,1:end-Delta[v]], label="rehat_$v")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("External input")
# end    
# # ----    


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

# # ----    
# # Plot the predicted internal inputs
# for v in 1:n
#     figure("Predicted internal input", figsize=(10,5))
#     PyPlot.plot(t, ri[v,:], label="ri_$v(t)")
#     PyPlot.plot(t[Delta[v]+1:end], rihat[v,1:end-Delta[v]], label="rihat_$v(t)")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Predicted internal input")
# end    
# # ----    


# # ----    
# # Plot the service rate, arrival rate, and addmitance rate
# for v in 1:n
#     figure("Service, arrival, and admittance rate - F_$v", figsize=(10,5))
#     PyPlot.plot(t, s[v,:]./dt, label="s_$v(t)")
#     PyPlot.plot(t, a[v,:]./dt, label="a_$v(t)")
#     PyPlot.plot(t, r[v,:]./dt, label="r_$v(t)")
#     xlabel("time (s)")
#     ylabel("packets per second")
#     legend(loc="lower right")
#     title("Service, arrival, and admittance rate - F_$v")
# end    
# # ----    

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

# # ----    
# # Plot average utility
# figure("Average efficiency ", figsize=(10,5))
# PyPlot.plot(t, Ue, label="Ue(t)")
# xlabel("time (s)")
# ylabel("packets per second")
# legend(loc="lower right")
# title("Average efficiency")

# figure("Average availability ", figsize=(10,5))
# PyPlot.plot(t, Ua, label="Ua(t)")
# xlabel("time (s)")
# ylabel("packets per second")
# legend(loc="lower right")
# title("Average availability")

# figure("Average utility function ", figsize=(10,5))
# PyPlot.plot(t, U, label="U(t)")
# xlabel("time (s)")
# ylabel("packets per second")
# legend(loc="lower right")
# title("Average utility function")
# # ----    

# ----
# Plot the efficiency
for v in 1:n
    figure("Efficiency", figsize=(10,5))
    PyPlot.plot(t, ue[v,:], label="e_$v(t)")
    xlabel("time (s)")
    ylabel("utility")
    legend(loc="lower right")
    title("Efficiency")
end    
# ----    
