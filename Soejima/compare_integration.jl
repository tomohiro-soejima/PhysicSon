include("./Find_average_N1.jl")
using .FindN1
include("./efficient_integral.jl")
using .integrateN1

x = [integrateN1.integrate_N1(10^i,5,10^-1) for i in 2 : 9]
y = [3.0/4*i*log(10) for i in 2 : 9]
println(x)
println(y)
