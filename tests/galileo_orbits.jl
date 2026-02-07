include("../orbits.jl")
include("../kepler.jl")
include("../reltools.jl")

using LinearAlgebra
using Plots
import .Orbits
import .RelTools

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
    GalileoHeight = 29599.8e+3
    GalileoEcc = 0.0
    GalileoInc = 56.0
    A_RAAN = 317.632
    B_RAAN = 77.632
    C_RAAN = 197.632
    R_Earth = 6.371e+6
    
    #=
    
        1 Galileo orbit = 67897.2231 seconds

    =#


    GSAT0218 = Orbits.new_orbit(GalileoHeight + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0) # A01
    GSAT0220 = Orbits.new_orbit(GalileoHeight + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B01
    GSAT0214 = Orbits.new_orbit(GalileoHeight + R_Earth, GalileoEcc, GalileoInc, C_RAAN, 0) # C01
    GSAT0226 = Orbits.new_orbit(GalileoHeight + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0) # A02
    #GSAT0221 = Orbits.new_orbit(GalileoHeight + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B02
    print(GSAT0218.T)
    #times = range(0.0, 67897.0, step=60) # 1 complete orbit in seconds, evaluated every minute

    plot(GSAT0218.x.(times), GSAT0218.y.(times), GSAT0218.z.(times)) # plotting GSAT0218
    #plot(times, GSAT0218.x.(times))
    #plot(GSAT0218.r.(times))
    #plot([GSAT0218.x.(times), GSAT0220.x.(times), GSAT0214.x.(times)], [GSAT0218.y.(times), GSAT0220.y.(times), GSAT0214.y.(times)], [GSAT0218.z.(times), GSAT0220.z.(times), GSAT0214.z.(times)]) # plotting GSAT0218
    

end

main()