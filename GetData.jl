#=

Simple script to extract randomly retrive one of the 1000 1h-data
files

Victor Millnert
Department of Automatic Control Theory
Lund University, Sweden

=#
module GetData

using JSON
using PyPlot
using Interpolations

using Distributions



"""

function getInputData(dt_fine::Float64, peak::Float64, plot_data::Bool)

It randomly selects one of the 1000 1-hour data files. The output is a
arrival-vector and a corresponding time-vector
    
"""
function getInputData(sim_time::Int64,
                      dt_fine::Float64,
                      peak::Float64,
                      stochastic::Bool,
                      plot_data::Bool)



    # randomly select one of the 100 input-files
    if sim_time == 1
        i = rand(1:1000)
        file = JSON.parsefile("./input_1h/input_$i.txt")
    elseif sim_time == 2
        i = rand(1:100)
        file = JSON.parsefile("./input_2h/input_$i.txt")        
    elseif sim_time == 6
        i = rand(1:100)
        file = JSON.parsefile("./input_6h/input_$i.txt")
    end

    # extract the "data"-fields, containing the number of packets
    data = file["data"];
    dt = file["meta"]["step"];
    t_start = file["meta"]["start"];
    t_end = file["meta"]["end"];
    # println("step size is $(dt) ")

    N = length(data)

    # convert the data-array from vector of {Any,1} to a vector of Float64
    input = Vector{Float64}(N)
    for i in 1:N
        if typeof(data[i][1]) == Float64
            input[i] = data[i][1];
        else
            input[i] = 0.0
        end
    end

    # convert the starting-time and stop-time to a time vector
    t = t_start:dt:t_end;
    if length(t) != N
        t = t_start:dt:t_end-dt;
    end

    # interpolate the input
    itp = interpolate(input, BSpline(Cubic(Flat())), OnCell());

    i_fine = 1:dt_fine/dt:N;

    N_fine = length(i_fine)
    t_fine = t_start:dt_fine:t_end

    N_min = min(length(i_fine), length(t_fine))
    t_fine = t_fine[1:N_min]
    i_fine = i_fine[1:N_min]
    
    # if length(i_fine) != length(t_fine)
    #     i_fine = 1:dt_fine/dt:N;
    #     t_fine = t_start:dt_fine:t_end-dt;
    # end

    # if mod(N_fine/10,10) != 0
    #     i_fine = i_fine[1:end-1]
    #     t_fine = t_fine[1:end-1]
    # end

    fine = zeros(length(i_fine))
    for i in 1:length(i_fine)
        fine[i] = itp[i_fine[i]]
    end
#    fine = itp[i_fine];


    # scale the input to make it a bit more changing!
    fine = fine - minimum(fine)/1.1
    input = input
    # scale the input such that the peak is at the specified peak
    # scale the input 
    scale = peak/maximum(fine)
    fine = fine*scale
    input = input*scale

    fine_rand = zeros(length(i_fine))
    # add stochasticity to the input
    if stochastic
        for i in 1:length(i_fine)
            fine_rand[i] = rand(Weibull(peak/100, fine[i]))
        end
    end
    
    if plot_data
        
        # ---------------
        # plot the output
        # ---------------
        t_date = map(x-> Dates.unix2datetime(x), t);
        t_fine_date = map(x-> Dates.unix2datetime(x), t_fine);
        
        close("all")
        # figure("Incoming octets", figsize=(15,10))
        # plot(t_date, input)

        if !stochastic
            figure("Interpolated data", figsize=(15,10))
            plot(t_fine, fine, label="interpolated")
            legend(loc="upper right")
        else
            figure("Comparison data", figsize=(15,10))
            plot(t_fine_date, fine_rand, label="stochastic")
            # plot(t_date, input, label="sampled")
            plot(t_fine_date, fine, label="interpolated")
            legend(loc="upper right")
        end
    end

    t_fine = t_fine-t_fine[1]
    if stochastic 
        return fine_rand, t_fine
    else
        return fine, t_fine
    end
    
end #end function

end #end module
