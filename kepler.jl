module Kepler

#=
    Kepler's Equation
=#
function K(t, a, ϵ, μ, E = π)
    M = sqrt(μ/(a^3)) * t
    return E - ϵ*sin(E) - M
end


#=
    Derivative of Kepler's Equation with respect to E
=#
function K_diff(t, a, ϵ, μ, E = π)
    M = sqrt(μ/(a^3)) * t
    return 1 - ϵ*cos(E)
end

#=
    This function takes in time as a parameter then uses Newton-Raphson method to solve Kepler's Equation.
=#
function E(t, a, ϵ, μ, E_prev = π, n = 0)
    if n == 5
        # Stop at five iterations
        return rad2deg(E_prev)
    end

    if abs(K(t, a, ϵ, μ, E_prev)) <= 0.001
        # Or stop when Kepler's equation is sufficiently close to zero
        return rad2deg(E_prev)
    else
        # Continue to the next iteration if neither stopping condition has been met
        E_next = E_prev - K(t, a, ϵ, μ, E_prev)/K_diff(t, a, ϵ, μ, E_prev)
        n += 1
        E(t, a, ϵ, μ, E_next, n)
    end
end

end