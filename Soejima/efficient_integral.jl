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
function integrate_N1_old(Length, n, v_min)
    v_max = 1 + v_min
    max_m = Int(floor(Length * (1 / v_min - 1 / v_max)))
    u_list = zeros(Float64, max_m + 2)
    N_list = zeros(Float64, max_m + 2)
    u_list[1] = v_min

    give_lists2!(u_list, N_list, max_m, v_min, n, Length)

    u_list[max_m + 2] = v_max
    N_list[max_m + 2] = 1

    p_here(x) = p(x, n, v_min)
    N_list .= p_here.(u_list) .* N_list

    #=
    println("F_list[max_m+1] = ", F_list[max_m+1])
    println("Np_list contains NaN:", any(isnan,Np_list))
    println("u_list contains NaN:", any(isnan,u_list))
    =#

    area = integrate(u_list, N_list)
    return area, u_list, N_list
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

function give_lists1!(u_list, N_list,max_m,v_min,n, Length)
    #not working yet
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
end

function give_lists2!(u_list,  N_list, max_m, v_min, n, Length)
    for m in reverse(0 : max_m)
        u_list[m + 1] = (1 / v_min - m / Length) ^ -1
        if m == max_m
            N_list[m + 1] = 1
        else
            N_list[m + 1] = N_list[m + 2] * (1 - F(u_list[m+2], n, v_min)) + 1
        end
    end
end

end #for module
