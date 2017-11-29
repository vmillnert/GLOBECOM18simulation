################################################
# ---------------------------------------------#
# ---------------------------------------------#
#          Plot the simulation results         #
# ---------------------------------------------# 
# ---------------------------------------------#
################################################

using PyPlot

# Plot metrics as a function of lambda
close("all")

# deadline_ratio vs. utility
figure("deadline_ratio vs. utility", figsize=(10,5))
PyPlot.scatter(deadline_ratio, u_mean, label="mean utility")
xlabel("deadline_ratio")
ylabel("Utility")
legend(loc="lower right")
title("deadline_ratio vs. utility")
# ----    

# ----
# deadline_ratio vs. discarded packets and overallocation
figure("deadline_ratio vs. discarded packets & overallocation", figsize=(10,5))
PyPlot.scatter(deadline_ratio, discarded*100, label="% discarded")
PyPlot.scatter(deadline_ratio, overallocation*100, label="% overallocation")
xlabel("deadline_ratio")
ylabel("percent")
legend(loc="lower right")
title("deadline_ratio vs. discarded packets & overallocation")
# ----    

# ----
# deadline_ratio vs. efficiency & availability
figure("deadline_ratio vs. efficiency", figsize=(10,5))
# PyPlot.plot(lam, U_mean, label="mean utility")
PyPlot.scatter(deadline_ratio, e_mean, label="efficiency")
PyPlot.scatter(deadline_ratio, a_mean, label="availability")
xlabel("deadline_ratio")
ylabel("percent")
legend(loc="lower right")
title("deadline_ratio vs. efficiency & availability")
# ----    
