using DifferentialEquations
using Plots
pyplot()

function N_body_Kuramoto!(du,u,p,t) #微分dy、関数u、パラメータp、変数t
    K = p[1]
    N = Int(p[2])
    omega = p[3:N+2]
    interaction = p[N+3:end]
    Kuramoto_interaction!(u,interaction)
    du .= omega .+ K.*interaction
end

function Kuramoto_interaction!(u,interaction)
    for i in 1:length(u)
        interaction = sum.(sin.([u[i]-u[j] for j in 1:length(u)]))
    end
end

u0 = [0,pi]
tspan = (0.0,100.0)
p = [0,2,2pi,2.2pi,0,0]
@time prob = ODEProblem(N_body_Kuramoto!,u0,tspan,p)
@time sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8,saveat=0.1)
