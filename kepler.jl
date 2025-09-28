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
        return E_prev * (180/π)
    end

    E_next = E_prev - K(t, a, ϵ, E_prev)/K_diff(t, a, ϵ, E_prev)
    relc = abs((E_next - E_prev)/E_prev) # the relative change from E_prev to E_next

    # stop the algorithm once successive iterations are within 2% of each other
    if relc <= 0.02
        return E_next * (180/π)
    else
        n += 1
        E(t, a, ϵ, E_next, n)
    end
    
end

end