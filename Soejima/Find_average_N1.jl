module FindN1

using LinearAlgebra
using QuadGK

export N1

"""
The probability distribution for car velocity, defined in paper. v_min is the minimum velocity and n is the order of F(v), also defined in paper.
"""
function p(v,n::Int,vmin)
    if vmin <= v <= 1+vmin
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

end #for module
