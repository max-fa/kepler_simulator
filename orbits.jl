module Orbits

export kepler_constructor

#=

This struct will contain the functions giving the
cartesian coordinates for an orbit

=#
struct Orbit
    x
    y
    z
end

#=

a = semimajor axis of orbit in meters
ϵ = eccentricity of orbit
i = angle of inclination of orbital plane with respect to the equatorial plane
Ω = longitude of the ascending node
ω = argument of the perigee

kepler_constructor() returns an Orbit{} struct

=#
function kepler_constructor(a::Real, ϵ::Real, i::Real, Ω::Real, ω::Real)

    r = φ::Real->a*(1-ϵ^2)/(1+ϵ*cosd(φ)) # r(φ), as derived from central force theory
    x0 = φ::Real->r(φ)*cosd(φ) # cartesian x
    y0 = φ::Real->r(φ)*sind(φ) # cartesian y

    # these are the x, y, and z components after applying rotation matrix
    x = φ::Real->(x0(φ) * (cosd(Ω)*cosd(ω) - sind(Ω)*sind(ω)) + y0(φ) * (sind(Ω)*cosd(i)*cosd(ω) + cosd(Ω)*sind(ω)))
    y = φ::Real->(x0(φ) * (-sind(Ω)*cosd(ω) - cosd(Ω)*cosd(i)*sind(ω)) + y0(φ) * (cosd(Ω)*cosd(i)*cosd(ω) - sind(Ω)*sind(ω)))
    z = φ::Real->(x0(φ) * sind(i)*sind(ω) - y0(φ) * sind(i)*cosd(ω))

    
    return Orbit(x, y, z)
end

end