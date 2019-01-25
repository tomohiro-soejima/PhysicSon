include("./Find_average_N1.jl")
using .FindN1
using Plots
using LaTeXStrings
pyplot()


N1_exp = zeros(Float64,6)
xaxis = [1,10,100,1000,10000,100000]


@time N1_exp[1] = N1(1,1,10^-2)[1]
@time N1_exp[2] = N1(10,1,10^-2)[1]
@time N1_exp[3] = N1(100,1,10^-2)[1]
@time N1_exp[4] = N1(1000,1,10^-2)[1]
@time N1_exp[5] = N1(10000,1,10^-2)[1]
@time N1_exp[6] = N1(100000,1,10^-2)[1]

plot(xaxis,N1_exp)
