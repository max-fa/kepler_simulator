include("../orbits.jl")
include("../kepler.jl")
include("../reltools.jl")

using cereal
using LinearAlgebra
using Plots
import .Orbits
import .RelTools

#=

Function for computing where in its orbit a satellite should be in
order to be on the past light cone of a specified user test location

=#
function emission_inverse(X_test, orbit)
    # [blank]
end

#=

Lunar Osculating Orbit Elements:
    semimajor axis = a = 384400e+3m
    eccentricity = ϵ = 0.05490
    inclination = i = 5.145°
    longitude of ascending node = Ω = 348°
    argument of the perigee = ω = 67.6°

=#


function main()
    
    # constellation-wide orbital parameters for Galileo satellites
    GalileoRadius = 29599.8e+3
    GalileoEcc = 0.0
    GalileoInc = 56.0
    A_RAAN = 317.632
    B_RAAN = 77.632
    C_RAAN = 197.632
    R_Earth = 6.371e+6
    locator = cereal.locatorselect(4, "FHC22")


    GSAT0218 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0) # A01
    GSAT0220 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B01
    GSAT0214 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, C_RAAN, 0) # C01
    GSAT0226 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0) # A02
    #GSAT0221 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B02
    
    # define test user location
    X_test = (time, 0, 0, R_Earth)

    X1 = emission_inverse(X_test, GSAT0218)
    X2 = emission_inverse(X_test, GSAT0220)
    X3 = emission_inverse(X_test, GSAT0214)
    X4 = emission_inverse(X_test, GSAT0226)

end

main()