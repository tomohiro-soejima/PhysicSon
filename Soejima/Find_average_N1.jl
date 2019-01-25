using LinearAlgebra
using NLsolve
using Roots
using Plots
using LaTeXStrings
using QuadGK
pyplot()

"""
The probability distribution for car velocity, defined in paper. v_min is the minimum velocity and n is the order of F(v), also defined in paper.
"""
function p(v,n::Int,vmin)
    if vmin< v < 1+vmin
        return n*(v-vmin)^(n-1)
    else
        return 0
    end
end


"""
Cumulative probability distribution of car velocity
"""
function F(v,n::Int,vmin)
    if vmin<v<1
        return (v-vmin)^n
    elseif v>1+vmin
        return 1
    else
        return 0
    end
end

function F̃(v_inv,n::Int,vmin)
    return F(1/v_inv,n,vmin)
end

Δt = 1
Length = 1

"""
Calculate the expected size of the first chunk by evaluating the sum and the product series explicitly
"""
function N(v,Length,n,v_min)
    N = 0
    curprod = 1
    for m in 1:(Int(round(Length/v))+1)
        N += curprod
        curprod *= (1-F̃(1/v - m/Length,n,v_min))
    end
    return N
end

"""
Calculate the expected value of <N_1> by integrating over N(v,L)*P(v)
"""
function N1(Length,n,v_min)
    N1 = quadgk(v->N(v,Length,n,v_min)*p(v,n,v_min),v_min,v_min+1)
    return N1
end

@time N1_1 = N1(1,1,10^-5)
@time N1_10 = N1(10,1,10^-5)
@time N1_100 = N1(100,1,10^-5)
@time N1_1000 = N1(1000,1,10^-5)

println(N1_1)
println(N1_10)
println(N1_100)
println(N1_1000)

@time N1_1_2 = N1(1,2,10^-5)
@time N1_10_2 = N1(10,2,10^-5)
@time N1_100_2 = N1(100,2,10^-5)
@time N1_1000_2 = N1(1000,2,10^-5)

println(N1_1_2)
println(N1_10_2)
println(N1_100_2)
println(N1_1000_2)

@time N1_1_2 = N1(1,20,10^-5)
@time N1_10_2 = N1(10,20,10^-5)
@time N1_100_2 = N1(100,20,10^-5)
@time N1_1000_2 = N1(1000,20,10^-5)

println(N1_1_20)
println(N1_10_20)
println(N1_100_20)
println(N1_1000_20)

println("end")

#=
f1(t) = N₁_prob(t,1)
f10(t) = N₁_prob(t,10)
f100(t) = N₁_prob(t,100)
f1000(t) = N₁_prob(t,1000)

x = 10^-5:0.04*10^-5:4*10^-5
y1 = f1.(x)
y10 = f10.(x)
y100 = f100.(x)
y1000 = f1000.(x)


plot(x,y1,size = (500,500),label = "L = 1",linewidth = 3)
plot!(x,y10,label="L = 10",linewidth = 3)
plot!(x,y100, label="L = 100",linewidth = 3)
plot!(x,y1000, label="L = 1000",linewidth = 3)
xlabel!("v",labelsize = 15)
ylabel!("N(v,L)",fontsize=15)
title!("N(v,L) vs v")
plot!(legendfontsize = 15,titlefontsize = 20,tickfontsize = 12,guidefontsize = 20)
savefig("N_v.png")
=#
