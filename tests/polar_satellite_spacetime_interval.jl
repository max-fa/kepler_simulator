include("../orbits.jl")
include("../kepler.jl")
include("../reltools.jl")

using cereal
using LinearAlgebra
using Plots
import .Orbits
import .RelTools

function main()
    
    # constellation-wide orbital parameters for Galileo satellites
    GalileoRadius = 29599.8e+3
    GalileoEcc = 0.0
    GalileoInc = 56.0
    A_RAAN = 317.632
    B_RAAN = 77.632
    C_RAAN = 197.632
    R_Earth = 6.371e+6
    c = 3e+8


    GSAT0218 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0) # A01
    GSAT0220 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B01
    GSAT0214 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, C_RAAN, 0) # C01
    GSAT0226 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0) # A02
    #GSAT0221 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B02
    
    #
    user_spherical = [R_Earth, 0, 0] # set in north pole in polar coordinates (r,θ,φ)
    user_cartesian = [0, 0, R_Earth] # north pole in cartesian coordinates (x,y,z)
    lightlikes = 0
    spacelikes = 0
    timelikes = 0
    for time in range(0.0, 67897.0, step=1)
        X = [time, GSAT0218.x(time), GSAT0218.y(time), GSAT0218.z(time)]
        user = [[time]; user_cartesian]
        S = X - user
        interval = -(c*S[1])^2 + S[2]^2 + S[3]^2 + S[4]^2
        #=print(interval)
        print("\n")=#
        if interval < 0
            timelikes += 1
        elseif interval > 0
            spacelikes += 1
        else
            lightlikes += 1
        end
    end

    print("Lightlikes: $(lightlikes)\n")
    print("Spacelikes: $(spacelikes)\n")
    print("Timelikes: $(timelikes)\n")


end

main()