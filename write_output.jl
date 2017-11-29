#=
Simple script to save down the output of the simulation to a text-file
Victor Millnert at the Department of Automatic Control
Theory, Lund University
=#


function save_long_sim_data(Delta_max::Int64, t, Rpath, r, scap, m, xi, U_mean, a)

    # how many data points do we need?
    data_points = 5e3 # 1000 data points
    dt = round(Int64,length(t)/1e3)
    
    i_start = round(Int64, 2*Delta_max)
    t = t[i_start:end]
    t = t - t[1]
    Rpath = Rpath[:,i_start:end]
    r = r[:,i_start:end]
    scap = scap[:,i_start:end]
    m = m[:,i_start:end]
    xi = xi[:,i_start:end]
    U_mean = U_mean[i_start:end]
    a = a[:,i_start:end]
    
    t = t./(60*60)
    folder_name = "data_long_sim_$(now())"
    run(`mkdir $(folder_name)`)

    f = open("$(folder_name)/data.txt", "w")
    write(f, "time  Rpath_1 Rpath_2 Rpath_3 r_1 r_2 r_3 scap_1 scap_2 scap_3 m_1 m_2 m_3 xi_1 xi_2 xi_3 U_mean a_1 a_2 a_3 \n")

    for i in 1:dt:length(t)
        write(f, "$(t[i]) $(Rpath[1,i]) $(Rpath[2,i]) $(Rpath[3,i]) $(r[1,i]) $(r[2,i]) $(r[3,i]) $(scap[1,i]) $(scap[2,i]) $(scap[3,i]) $(m[1,i]) $(m[2,i]) $(m[3,i]) $(xi[1,i]) $(xi[2,i]) $(xi[3,i]) $(U_mean[i]) $(a[1,i]) $(a[2,i]) $(a[3,i]) \n")
    end
    close(f)    

end #end function

function save_lambda_sim_data(lambda,
                              u_mean,
                              e_mean,
                              a_mean,
                              discarded,
                              overallocation)

    N = length(lambda)
    
    folder_name = "data_lambda_sim_$(now())"
    run(`mkdir $(folder_name)`)

    # write the path-wise input rate
    f = open("$(folder_name)/u_mean.txt", "w")
    for i in 1:N
        write(f, "$(lambda[i]) \t $(u_mean[i]) \t \n")
    end
    close(f)    

    # write the node-wise input rate
    f = open("$(folder_name)/e_a.txt", "w")
    for i in 1:N
        write(f, "$(lambda[i]) \t $(e_mean[i]) \t $(a_mean[i]) \t \n")
    end
    close(f)    

    # write the service capacity
    f = open("$(folder_name)/discarded_overallocation.txt", "w")
    for i in 1:N
        write(f, "$(lambda[i]) \t $(discarded[i]*100) \t $(overallocation[i]*100) \t \n")
    end
    close(f)    
end #end function

function save_deadline_ratio_sim_data(deadline_ratio, u_mean, e_mean, a_mean, discarded, overallocation,
                                      u_var, e_var, a_var, discarded_var, overallocation_var)

    N = length(deadline_ratio)
    
    folder_name = "data_deadline_ratio_sim_$(now())"
    run(`mkdir $(folder_name)`)

    # our data
    f = open("$(folder_name)/data.txt", "w")
    write(f, "ratio u_mean u_var e_mean e_var a_mean a_var discarded discarded_var  overallocation overallocation_var  \n")

    #     # write the path-wise input rate
    # f = open("$(folder_name)/u_mean.txt", "w")
    for i in 1:N
        write(f, "$(deadline_ratio[i]) $(u_mean[i]) $(u_var[i]) $(e_mean[i]) $(e_var[i]) $(a_mean[i]) $(a_var[i]) $(discarded[i]) $(discarded_var[i])  $(overallocation[i]) $(overallocation_var[i]) \n")
    end
    close(f)    

    # # write the node-wise input rate
    # f = open("$(folder_name)/e_a.txt", "w")
    # for i in 1:N
    #     write(f, "$(deadline_ratio[i]) \t $(e_mean[i]) \t $(a_mean[i]) \t \n")
    # end
    # close(f)    

    # # write the service capacity
    # f = open("$(folder_name)/discarded_overallocation.txt", "w")
    # for i in 1:N
    #     write(f, "$(deadline_ratio[i]) \t $(discarded[i]*100) \t $(overallocation[i]*100) \t \n")
    # end
    # close(f)    
end #end function

# function save_comparison_sim_data(u_mean, e_mean, a_mean, discarded, overallocation,
#                                   u_mean_DAS, e_mean_DAS, a_mean_DAS, discarded_DAS, overallocation_DAS,
#                                   u_mean_DAS_AC, e_mean_DAS_AC, a_mean_DAS_AC, discarded_DAS_AC, overallocation_DAS_AC,
#                                   u_mean_DOA, e_mean_DOA, a_mean_DOA, discarded_DOA, overallocation_DOA,
#                                   u_mean_DOA_AC, e_mean_DOA_AC, a_mean_DOA_AC, discarded_DOA_AC, overallocation_DOA_AC)
function save_comparison_sim_data(u_mean, e_mean, a_mean, discarded, overallocation,
                         u_var, e_var, a_var, discarded_var, overallocation_var,
                         u_mean_DAS, e_mean_DAS, a_mean_DAS, discarded_DAS, overallocation_DAS,
                         u_var_DAS, e_var_DAS, a_var_DAS, discarded_var_DAS, overallocation_var_DAS,
                         u_mean_DAS_AC, e_mean_DAS_AC, a_mean_DAS_AC, discarded_DAS_AC, overallocation_DAS_AC,
                         u_var_DAS_AC, e_var_DAS_AC, a_var_DAS_AC, discarded_var_DAS_AC, overallocation_var_DAS_AC,
                         u_mean_DOA, e_mean_DOA, a_mean_DOA, discarded_DOA, overallocation_DOA,
                         u_var_DOA, e_var_DOA, a_var_DOA, discarded_var_DOA, overallocation_var_DOA,
                         u_mean_DOA_AC, e_mean_DOA_AC, a_mean_DOA_AC, discarded_DOA_AC, overallocation_DOA_AC,
                         u_var_DOA_AC, e_var_DOA_AC, a_var_DOA_AC, discarded_var_DOA_AC, overallocation_var_DOA_AC)


    
    folder_name = "data_comparison_sim_$(now())"
    run(`mkdir $(folder_name)`)

    # our data
    f = open("$(folder_name)/data.txt", "w")
    write(f, "sim u_mean u_var e_mean e_var a_mean a_var discarded discarded_var  overallocation overallocation_var  \n")

    write(f, "1  $(u_mean) $(u_var) $(e_mean) $(e_var) $(a_mean) $(a_var) $(discarded) $(discarded_var)  $(overallocation) $(overallocation_var) \n")

    # DAS
    write(f, "2  $(u_mean_DAS) $(u_var_DAS) $(e_mean_DAS) $(e_var_DAS) $(a_mean_DAS) $(a_var_DAS) $(discarded_DAS) $(discarded_var_DAS)  $(overallocation_DAS) $(overallocation_var_DAS) \n")

    # write(f, "2  $(u_mean_DAS)  $(e_mean_DAS)  $(a_mean_DAS)  $(discarded_DAS)  $(overallocation_DAS)  \n")

    # DOA
    write(f, "3  $(u_mean_DOA) $(u_var_DOA) $(e_mean_DOA) $(e_var_DOA) $(a_mean_DOA) $(a_var_DOA) $(discarded_DOA) $(discarded_var_DOA)  $(overallocation_DOA) $(overallocation_var_DOA) \n")
    # write(f, "3  $(u_mean_DOA)  $(e_mean_DOA)  $(a_mean_DOA)  $(discarded_DOA)  $(overallocation_DOA)  \n")

    # DAS + AC
    write(f, "4  $(u_mean_DAS_AC) $(u_var_DAS_AC) $(e_mean_DAS_AC) $(e_var_DAS_AC) $(a_mean_DAS_AC) $(a_var_DAS_AC) $(discarded_DAS_AC) $(discarded_var_DAS_AC)  $(overallocation_DAS_AC) $(overallocation_var_DAS_AC) \n")
    # write(f, "4  $(u_mean_DAS_AC)  $(e_mean_DAS_AC)  $(a_mean_DAS_AC)  $(discarded_DAS_AC)  $(overallocation_DAS_AC)  \n")
    
    # DOA + AC
    write(f, "5  $(u_mean_DOA_AC) $(u_var_DOA_AC) $(e_mean_DOA_AC) $(e_var_DOA_AC) $(a_mean_DOA_AC) $(a_var_DOA_AC) $(discarded_DOA_AC) $(discarded_var_DOA_AC)  $(overallocation_DOA_AC) $(overallocation_var_DOA_AC) \n")
    # write(f, "5  $(u_mean_DOA_AC)  $(e_mean_DOA_AC)  $(a_mean_DOA_AC)  $(discarded_DOA_AC)  $(overallocation_DOA_AC)  \n")

    close(f)    

end #end function

