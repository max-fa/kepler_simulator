module Orbits

include("kepler.jl")
import .Kepler

μ = 3.99e+14

#=

This struct contains the classical Keplerian elements that define a unique geocentric orbit

=#
struct Elements
    a::Real
    ϵ::Real
    i::Real
    Ω::Real
    ω::Real
end

#=

This struct will contains the functions giving the polar or cartesian coordinates for an orbit,
as well as the Keplerian elements that define the specific orbit

=#
struct Orbit
    elements::Elements
    T::Real
    x::Function
    y::Function
    z::Function
    r::Function
    φ::Function
end


# returns cos(φ) at a given time, for a given orbit
function cosφ(t, a, ϵ)
    E = Kepler.E(t, a, ϵ)
    return (ϵ - cosd(E))/(ϵ*cosd(E) - 1)
end

# returns the orbital radius at a given time, for a given orbit
function r(t, a, ϵ)
    return a*(1 - ϵ^2)/(1 + ϵ*cosφ(t, a, ϵ))
end

# returns the orbital plane polar angle at a given time
function φ_t(t, a, ϵ)
    E = Kepler.E(t, a, ϵ)
    φ_raw = acosd((ϵ - cosd(E))/(ϵ*cosd(E) - 1)) # φ calculated from eccentric anomaly, always in [0,180]

    # cos(φ), I'm calculating directly from E instead of using the
    # cosφ() function so that I can avoid running Newton-Raphson again
    c = (ϵ - cosd(E))/(ϵ*cosd(E) - 1)

    if E in range(0,180)
        # if E is in the top hemisphere of the orbit, then φ will
        # be left in the default, upper hemisphere given by acosd
        return φ_raw
    else
        # E is in the lower hemisphere and φ is in the uppper.
        # Must compute the quadrant that φ is in in order to
        # properly shift it into the bottom hemisphere
        if c > 0
            # φ is in Quadrant I, shift to IV
            return 360 - φ_raw
        elseif c < 0
            # φ is in Quadrant II, so shift it to III
            return 180 + φ_raw
        else
            # φ is 90°
            return 270
        end
    end

end

# returns cartesian x coordinate at a given time, for a given orbit
function x_t(t, a, ϵ)
    E = Kepler.E(t, a, ϵ)
    return a * (cosd(E) - ϵ)
end

# returns cartesian y coordinate at a given time, for a given orbit
function y_t(t, a, ϵ)
    E = Kepler.E(t, a, ϵ)
    return a * sqrt(1 - ϵ^2) * sind(E)
end

#=

a = semimajor axis of orbit in meters
ϵ = eccentricity of orbit
i = angle of inclination of orbital plane with respect to the equatorial plane
Ω = longitude of the ascending node
ω = argument of the perigee

new_orbit() returns an Orbit struct

=#
function new_orbit(a::Real, ϵ::Real, i::Real, Ω::Real, ω::Real)
    
    elements = Elements(a, ϵ, i, Ω, ω)
    T = 2*π*sqrt((a^3)/μ)

    orbit_r = t->r(t, a, ϵ)
    orbit_φ = t->φ_t(t, a, ϵ)
    x0 = t->x_t(t, a, ϵ) # cartesian x
    y0 = t->y_t(t, a, ϵ) # cartesian y

    # these are the x, y, and z components after applying rotation matrix
    orbit_x = t->(x0(t) * (cosd(Ω)*cosd(ω) - sind(Ω)*sind(ω)) + y0(t) * (sind(Ω)*cosd(i)*cosd(ω) + cosd(Ω)*sind(ω)))
    orbit_y = t->(x0(t) * (-sind(Ω)*cosd(ω) - cosd(Ω)*cosd(i)*sind(ω)) + y0(t) * (cosd(Ω)*cosd(i)*cosd(ω) - sind(Ω)*sind(ω)))
    orbit_z = t->(x0(t) * sind(i)*sind(ω) - y0(t) * sind(i)*cosd(ω))
    

    
    return Orbit(elements, T, orbit_x, orbit_y, orbit_z, orbit_r, orbit_φ)
end

end