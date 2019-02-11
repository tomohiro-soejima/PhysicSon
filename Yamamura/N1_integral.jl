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
    p_here(x) = p(x, n, v_min)

    u_old = (1 / v_min - max_m / Length) ^ -1
    N_old = 1.0
    p_old = p_here(v_max)
    p_new = p_here(u_old)
    area = (v_max - u_old) * (p_new + p_old*N_old)/2
    p_old = p_new
    for m in reverse(1 : max_m)
        N_new = N_old * (1 - F(u_old, n, v_min)) + 1
        u_new = (1 / v_min - (m - 1) / Length) ^ -1
        p_new = p_here(u_new)
        area += (u_old - u_new) * (p_new * N_new + p_old * N_old)/2
        N_old = N_new
        u_old = u_new
        p_old = p_new
    end

    return area
end
end #for module
