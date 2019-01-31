include("./Find_average_N1.jl")
using .FindN1
using Plots
using LaTeXStrings
pyplot()

function calculateN1exp(number_of_data_points, order_of_Fv,vmin)
    N1_exp = zeros(Float64,number_of_data_points)
    xaxis = [10^(index-1) for index in 1:number_of_data_points]

    for (index, x_value) in enumerate(xaxis)
           @time N1_exp[index] = N1(x_value,order_of_Fv,vmin)[1]
    end

    println(N1_exp)

    return (xaxis,N1_exp)
end

N1_exp_10th_1 = calculateN1exp(7,10,10^-2)
#=
This was the result for N1_exp_10th_1 = calculateN1exp(7,10,10^-2)
0.001730 seconds (2.07 k allocations: 42.656 KiB)
0.060569 seconds (15.63 k allocations: 321.688 KiB)
0.065535 seconds (6.78 k allocations: 140.656KiB)
0.260737 seconds (6.68 k allocations: 138.781KiB)
5.601989 seconds (37.69 k allocations: 768.531 KiB)
321.008353 seconds (257.65 k allocations: 5.165MiB)
172349.480222 seconds (1.54 M allocations: 30.354 MiB, 0.00% gc time)
[1.09626, 1.38164, 3.2852, 6.17206, 12.8927, 57.2552, 477.982]
=#
