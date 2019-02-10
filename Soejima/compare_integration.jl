include("./Find_average_N1.jl")
using .FindN1
include("./efficient_integral.jl")
using .integrateN1

#FindN1.N1(1,1,10^-1)
#try1 = @time FindN1.N1(10,2,10^-3)
integrateN1.integrate_N1(1,1,10^-1)
integrateN1.integrate_N1_old(1,1,10^-1)
try1 = @time integrateN1.integrate_N1(1000000,2,10^-3)
try2, u_list, N_list = @time integrateN1.integrate_N1_old(1000000,2,10^-3)


println(try1)
println(try2)
