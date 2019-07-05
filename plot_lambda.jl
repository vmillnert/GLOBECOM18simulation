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

# lambda vs. utility
figure("Lambda vs. utility", figsize=(10,5))
PyPlot.scatter(lambda, u_mean, label="mean utility")
xlabel("lambda")
ylabel("Utility")
legend(loc="lower right")
title("Lambda vs. utility")
# ----    

# ----
# lambda vs. discarded packets and overallocation
figure("Lambda vs. discarded packets & overallocation", figsize=(10,5))
PyPlot.scatter(lambda, discarded*100, label="% discarded")
PyPlot.scatter(lambda, overallocation*100, label="% overallocation")
xlabel("lambda")
ylabel("percent")
legend(loc="lower right")
title("Lambda vs. discarded packets & overallocation")
# ----    

# ----
# lambda vs. efficiency & availability
figure("Lambda vs. efficiency", figsize=(10,5))
# PyPlot.plot(lam, U_mean, label="mean utility")
PyPlot.scatter(lambda, e_mean, label="efficiency")
PyPlot.scatter(lambda, a_mean, label="availability")
xlabel("lambda")
ylabel("percent")
legend(loc="lower right")
title("Lambda vs. efficiency & availability")
# ----    
