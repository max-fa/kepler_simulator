module Orbits

using DoubleFloats
include("kepler.jl")
import .Kepler
m_E = Double64(5.972e+24)
G = Double64(6.674e-11)
μ = G * m_E
c = Double64(2.998e+8)

#=

This struct contains the classical Keplerian elements that define a unique geocentric orbit

=#
struct Elements
    a::Double64
    ϵ::Double64
    i::Double64
    Ω::Double64
    ω::Double64
end

#=

This struct will contains the functions giving the polar or cartesian coordinates for an orbit,
as well as the Keplerian elements that define the specific orbit

=#
struct Orbit
    elements::Elements
    T::Double64
    x::Function
    y::Function
    z::Function
    r::Function
    φ::Function
    v::Function
    τ::Function
    name::String
    offset::Double64
end


# returns cos(φ) at a given time, for a given orbit
function cosφ(t, a, ϵ, offset)
    E = Kepler.E(t, a, ϵ, μ, offset)
    return (ϵ - cosd(E))/(ϵ*cosd(E) - Double64(1))
end

# returns the orbital radius at a given time, for a given orbit
function r(t, a, ϵ, offset)
    coeff = (Double64(1) - ϵ^Double64(2))/(Double64(1) + ϵ*cosφ(t, a, ϵ, offset))
    return a*coeff
end

function normalize(φ_raw, E, t, T)
    n = floor(E / Double64(180))
    
    return (((Double64(-1))^n)*(φ_raw) + (n + Double64(0.5)*(Double64(1) - (Double64(-1))^n))*Double64(180) - (Double64(360)*t/T)) * (Double64(180)/Double64(π))
end

# returns the orbital plane polar angle at a given time
function φ_t(t, a, ϵ, offset)
    T = Double64(2.0)*Double64(π)*sqrt((a^Double64(3))/μ)
    E = Kepler.E(t, a, ϵ, μ, offset)
    #return E
    φ_raw = acosd((ϵ - cosd(E))/(ϵ*cosd(E) - Double64(1))) # φ calculated from eccentric anomaly, always in [0,180]
    return φ_raw
    #return normalize(φ_raw, E, t, T)
end

# returns cartesian x coordinate at a given time, for a given orbit
function x_t(t, a, ϵ, offset)
    E = Kepler.E(t, a, ϵ, μ, offset)
    return a * (cosd(E) - ϵ)
end

# returns cartesian y coordinate at a given time, for a given orbit
function y_t(t, a, ϵ, offset)
    E = Kepler.E(t, a, ϵ, μ, offset)
    return a * sqrt(Double64(1) - ϵ^Double64(2)) * sind(E)
end

function v(t, a, ϵ, offset)
    _r = r(t, a, ϵ, offset)
    return sqrt(G*m_E*(Double64(2)/_r - Double64(1)/a))
end

function τ(t, a, ϵ, offset)
    γ = sqrt(Double64(1) - (v(t, a, ϵ, offset)^Double64(2))/(c^Double64(2)))
    #=print(γ)
    print("\n")
    print(t*γ)
    print("\n")
    print(t)
    print("\n\n")=#
    return γ*t
end

#=

a = semimajor axis of orbit in meters
ϵ = eccentricity of orbit
i = angle of inclination of orbital plane with respect to the equatorial plane
Ω = longitude of the ascending node
ω = argument of the perigee

new_orbit() returns an Orbit struct

=#
function new_orbit(a::Double64, ϵ::Double64, i::Double64, Ω::Double64, ω::Double64, name::String, offset::Double64=Double64(0))
    
    elements = Elements(a, ϵ, i, Ω, ω)
    T = Double64(2.0)*Double64(π)*sqrt((a^Double64(3))/μ)

    orbit_r = t->r(t, a, ϵ, offset)
    orbit_φ = t->φ_t(t, a, ϵ, offset)
    x0 = t->x_t(t, a, ϵ, offset) # cartesian x
    y0 = t->y_t(t, a, ϵ, offset) # cartesian y

    # these are the x, y, and z components after applying rotation matrix
    orbit_x = t->(x0(t) * (cosd(Ω)*cosd(ω) - sind(Ω)*sind(ω)) + y0(t) * (sind(Ω)*cosd(i)*cosd(ω) + cosd(Ω)*sind(ω)))
    orbit_y = t->(x0(t) * (-sind(Ω)*cosd(ω) - cosd(Ω)*cosd(i)*sind(ω)) + y0(t) * (cosd(Ω)*cosd(i)*cosd(ω) - sind(Ω)*sind(ω)))
    orbit_z = t->(x0(t) * sind(i)*sind(ω) - y0(t) * sind(i)*cosd(ω))

    orbit_v = t->v(t, a, ϵ, offset)
    orbit_τ = t->τ(t, a, ϵ, offset)

    
    return Orbit(elements, T, orbit_x, orbit_y, orbit_z, orbit_r, orbit_φ, orbit_v, orbit_τ, name, offset)
end

end