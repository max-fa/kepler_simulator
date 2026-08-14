include("../orbits.jl")
include("../cereal_inv.jl")
include("../reltools.jl")


using DoubleFloats
using Infiltrator
using cereal
using LinearAlgebra
using Plots
using .CerealInv
import .Orbits
import .RelTools

#pyplot()

function main()
    
    # constellation-wide orbital parameters for Galileo satellites
    GalileoRadius = Double64(29599.8e+3)
    GalileoEcc = Double64(0.0)
    GalileoInc = Double64(56.0)
    A_RAAN = Double64(317.632)
    B_RAAN = Double64(77.632)
    C_RAAN = Double64(197.632)
    R_Earth = Double64(6.371e+6)
    c = Double64(3e+8)
    A02_offset = Double64(3600)


    GSAT0218 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, Double64(0), "GSAT0218") # A01
    GSAT0220 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, Double64(0), "GSAT0220") # B01
    GSAT0214 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, C_RAAN, Double64(0), "GSAT0214") # C01
    GSAT0226 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, Double64(0), "GSAT0226", A02_offset) # A02
    #GSAT0221 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B02

    satellites = (GSAT0218, GSAT0220, GSAT0214, GSAT0226)
    
    reception_time = rand(range(Double64(0.0), Double64(86400), step=Double64(1)))
    #user_spherical = [R_Earth, 0, 0] # north pole in polar coordinates (r,θ,φ)
    user_cartesian = [reception_time, Double64(0), Double64(0), R_Earth] # north pole in cartesian coordinates (x,y,z)

    print("Reception time is: $(reception_time)\n\n")

    X_A01 = CerealInv.get_emission_vector(user_cartesian, GSAT0218)
    X_B01 = CerealInv.get_emission_vector(user_cartesian, GSAT0220)
    X_C01 = CerealInv.get_emission_vector(user_cartesian, GSAT0214)
    X_A02 = CerealInv.get_emission_vector(user_cartesian, GSAT0226)

    A01_emission_interval = RelTools.get_interval(X_A01 - user_cartesian)
    B01_emission_interval = RelTools.get_interval(X_B01 - user_cartesian)
    C01_emission_interval = RelTools.get_interval(X_C01 - user_cartesian)
    A02_emission_interval = RelTools.get_interval(X_A02 - user_cartesian)

    print("A01 emission time is: $(X_A01[1])\n\n")
    print("B01 emission time is: $(X_B01[1])\n\n")
    print("C01 emission time is: $(X_C01[1])\n\n")
    print("A02 emission time is: $(X_A02[1])\n\n")

    print("A01 separation spacetime interval at emission time is: $(A01_emission_interval)\n\n")
    print("B01 separation spacetime interval at emission time is: $(B01_emission_interval)\n\n")
    print("C01 separation spacetime interval at emission time is: $(C01_emission_interval)\n\n")
    print("A02 separation spacetime interval at emission time is: $(A02_emission_interval)\n\n")

    X = [X_A01 X_B01 X_C01 X_A02]

    locator = cereal.locatorselect(4, "FHC22")

    Xc = locator(X)

    print(Xc)

end

function count_positives(intervals)
    return sum(x->(x>0), intervals)
end

# returns an array of indices locating the right brackets of crossings in an intervals set
function find_brackets(intervals)
    i = 0
    brackets = []
    for interval in intervals
        i += 1
        if i != length(intervals)
            if intervals[i+1] * interval < 0
                push!(brackets, (i, i+1))
            end
        end
    end

    return brackets
end

main()