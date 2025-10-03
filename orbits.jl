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
    v::Function
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

function normalize(φ_raw, E, t, T)
    n = floor(E / 180)
    
    return (((-1)^n)*(φ_raw) + (n + 0.5*(1 - (-1)^n))*(180) - (360*t/T)) * (180/π)
end

# returns the orbital plane polar angle at a given time
function φ_t(t, a, ϵ)
    T = 2*π*sqrt((a^3)/μ)
    E = Kepler.E(t, a, ϵ)
    #return E
    φ_raw = acosd((ϵ - cosd(E))/(ϵ*cosd(E) - 1)) # φ calculated from eccentric anomaly, always in [0,180]
    return normalize(φ_raw, E, t, T)
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

function v(t, a, ϵ)
    _r = r(t, a, ϵ)
    return sqrt(μ*(2/_r - 1/a))
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

    orbit_v = t->v(t, a, ϵ)
    

    
    return Orbit(elements, T, orbit_x, orbit_y, orbit_z, orbit_r, orbit_φ, orbit_v)
end

end