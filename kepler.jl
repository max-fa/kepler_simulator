module Kepler

using ForwardDiff

E_counter = 0
#=
    This function is Kepler's Equation.
    E(t) will be finding the roots of this function by varying E.
=#
function K(E, ϵ, a, t)
    μ = 3.99e+14 # assuming geocentric orbit
    M = sqrt(μ/(a^3)) * t
    return E - ϵ*sin(E) - M
end

function K_diff(E, ϵ, a, t)
    μ = 3.99e+14
    M = sqrt(μ/(a^3)) * t
    return 1 - ϵ*cos(E)
end

#=
    This function takes in time as a parameter then uses Newton-Raphson method to solve Kepler's Equation.
    E0 is the guess for E to be used in the Newton-Raphson algorithm.
    This gives the eccentric anomaly, E, that can then be used to find the radius and true anomaly.
=#
function E(E_prev, ϵ, a, t)
    
    global E_counter += 1
    E_next = E_prev - K(E_prev, ϵ, a, t)/K_diff(E_prev, ϵ, a, t)
    relc = abs((E_next - E_prev)/E_prev) # the relative change from E_prev to E_next

    if E_counter == 5
        return E_next
    end

    relc <= 0.02 ? (return E_next) : E(E_next, ϵ, a, t) # stop the algorithm once successive iterations are within 2% of each other
    
end

end