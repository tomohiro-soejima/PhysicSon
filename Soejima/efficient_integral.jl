module integrateN1

using LinearAlgebra

export integrate_N1

"""
The probability distribution for car velocity, defined in paper. v_min is the minimum velocity and n is the order of F(v), also defined in paper.
"""
function p(v,n::Int,vmin)
    if vmin <= v <= 1+vmin
        return n*(v-vmin)^(n-1)
    else
        return zero(v)
    end
end


"""
Cumulative probability distribution of car velocity
"""
function F(v,n::Int,vmin)
    if vmin<=v<=1+vmin
        return (v-vmin)^n
    elseif v>=1+vmin
        return oneunit(v)
    else
        return zero(v)
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
function integrate_N1(Length, n, v_min)
    v_max = 1 + v_min
    max_m = Int(floor(Length * (1 / v_min - 1 / v_max)))
    u_list = zeros(Float64, max_m + 1)
    F_list = zeros(Float64, max_m + 1)
    N_list = zeros(Float64, max_m + 1)
    Np_list = zeros(Float64, max_m + 1)
    u_list[1] = v_min
    F_list[1] = 1
    cur_prod = 1
    N = 0
    N += cur_prod
    for m in 1 : max_m
        u_list[m + 1] = (1 / v_min - m / Length) ^ -1
        F_list[m + 1] =  1 - F(u_list[m + 1], n, v_min)
        cur_prod *= F_list[m + 1]
        N += cur_prod
    end

    N_list[1] = N

    for index in 1 : max_m
        N_list[index + 1] = (N_list[index] - 1.0) / F_list[index + 1]
    end

    p_here(x) = p(x, n, v_min)

    Np_list = p_here.(u_list) .* N_list

    println("F_list[max_m+1] = ", F_list[max_m+1])
    println("Np_list contains NaN:", any(isnan,Np_list))
    println("u_list contains NaN:", any(isnan,u_list))

    area = integrate(u_list, Np_list)
    return area, u_list, Np_list, F_list
end

function integrate(u_list,Np_list)
    M =length(u_list)
    area = 0
    for index in 1:M-1
        du = u_list[index+1] - u_list[index]
        area += du * (Np_list[index + 1] + Np_list[index]) / 2
    end
    return area
end


end #for module
