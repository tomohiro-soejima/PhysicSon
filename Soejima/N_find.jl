using LinearAlgebra
using NLsolve
using Roots
using Plots
pyplot()

function f₀(v)
    if 10^-5< v < 1+10^-5
        return 1
    else
        return 0
    end
end

function L(v)
    if 1+10^-5>v>10^-5
        return (v-10^-5)
    elseif v>1+10^-5
        return 1
    else
        return 0
    end
end

function L̃(v)
    return L(1/v)
end

function return_v()
    a = rand()
    find_zero(v->L(v)-a,10^-5)
end

Δt = 1
Length = 1

function chunk_average(Length,N)
    N₁ = zeros(Int64,N)
    for i in 1:N
        N₁[i] = chunk_size(Length)
    end
    max = maximum(N₁)
    average = sum(N₁)/N
    variance = sqrt(sum(N₁.^2))/N
    return average,variance,max
end

function chunk_size(Length)
    v₁ = return_v()
    return N₁_seq(v₁,Length)
end

function N₁_seq(v₁,Length)
    n = 1
    while true
        vᵢ = return_v()
        if vᵢ>L̃(1/v₁ - n/Length)
            n += 1
        else
            return n
        end
    end
end

function N₁_prob(v,Length)
    N = 0
    curprod = 1
    for m in 1:(Int(round(Length/v))+1)
        N += curprod
        curprod *= (1-L̃(1/v - m/Length))
    end
    return N
end

#=
x = 1:100
f(t) = N₁_prob(1*10^-5,t)
y = f.(x)
plot(x,y,size =(200,200))=#
#plot(x,y./sqrt.(x),size = (200,200))
##plot(x.^(1/8),y,size=(200,200))
#plot(log10.(x),log10.(y),size =(200,200))

x = 10^-5:0.04*10^-5:4*10^-5
f1(t) = N₁_prob(t,1)
f10(t) = N₁_prob(t,10)
f100(t) = N₁_prob(t,100)
f1000(t) = N₁_prob(t,1000)
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

plot(x,1 ./ y1,size = (500,500),label = "L = 1",linewidth = 3)
plot!(x,1 ./ y10,label="L = 10",linewidth = 3)
plot!(x,1 ./ y100, label="L = 100",linewidth = 3)
plot!(x,1 ./ y1000, label="L = 1000",linewidth = 3)
xlabel!("v",labelsize = 15)
ylabel!("X(v)^-1",fontsize=15)
title!("X(v)^-1 vs v")
plot!(legendfontsize = 15,titlefontsize = 20,tickfontsize = 12,guidefontsize = 20)
savefig("N_v_inv.png")

x = 10^-5:0.01*0.01:0.01
f1(t) = N₁_prob(t,1)
f10(t) = N₁_prob(t,10)
f100(t) = N₁_prob(t,100)
f1000(t) = N₁_prob(t,1000)
y1 = f1.(x)
y10 = f10.(x)
y100 = f100.(x)
y1000 = f1000.(x)
plot(x,1 ./ y1,size = (500,500),label = "L = 1",linewidth = 3)
plot!(x,1 ./ y10,label="L = 10",linewidth = 3)
plot!(x,1 ./ y100, label="L = 100",linewidth = 3)
plot!(x,1 ./ y1000, label="L = 1000",linewidth = 3)
xlabel!("v",labelsize = 15)
ylabel!("X(v)^-1",fontsize=15)
title!("X(v)^-1 vs v")
plot!(legendfontsize = 15,titlefontsize = 20,tickfontsize = 12,guidefontsize = 20)
savefig("X_v_inv_2.png")


x = 10^-5:0.5*10^-6:10^-5+10*10^-6
f1(t) = N₁_prob(t,1)
f10(t) = N₁_prob(t,10)
f100(t) = N₁_prob(t,100)
f1000(t) = N₁_prob(t,1000)
y1 = f1.(x)
y10 = f10.(x)
y100 = f100.(x)
y1000 = f1000.(x)
plot(x,1 ./ y1,size = (500,500),label = "L = 1",linewidth = 3)
plot!(x,1 ./ y10,label="L = 10",linewidth = 3)
plot!(x,1 ./ y100, label="L = 100",linewidth = 3)
plot!(x,1 ./ y1000, label="L = 1000",linewidth = 3)
xlabel!("v",labelsize = 15)
ylabel!("X(v)^-1",fontsize=15)
title!("X(v)^-1 vs v")
plot!(legendfontsize = 15,titlefontsize = 20,tickfontsize = 12,guidefontsize = 20)
savefig("X_v_inv_3.png")
