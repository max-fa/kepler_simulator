module Kepler

μ = 3.99e+14

#=
    Kepler's Equation
=#
function K(t, a, ϵ, E = π)
    M = sqrt(μ/(a^3)) * t
    return E - ϵ*sin(E) - M
end


#=
    Derivative of Kepler's Equation with respect to E
=#
function K_diff(t, a, ϵ, E = π)
    M = sqrt(μ/(a^3)) * t
    return 1 - ϵ*cos(E)
end

#=
    This function takes in time as a parameter then uses Newton-Raphson method to solve Kepler's Equation.
=#
function E(t, a, ϵ, E_prev = π, n = 0)
    if n == 5
        #=print("Imprecision")
        print("\n")=#
        return rad2deg(E_prev)
    end

    if abs(K(t, a, ϵ, E_prev)) <= 0.001
        #=print("Precision")
        print("\n")=#
        return rad2deg(E_prev)
    else
        E_next = E_prev - K(t, a, ϵ, E_prev)/K_diff(t, a, ϵ, E_prev)
        n += 1
        E(t, a, ϵ, E_next, n)
    end
end

end