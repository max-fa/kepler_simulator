include("../orbits.jl")
include("../cereal_inv.jl")
include("../reltools.jl")


using Infiltrator
using cereal
using LinearAlgebra
using Plots
using .CerealInv
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
    A02_offset = 3600


    GSAT0218 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0, "GSAT0218") # A01
    GSAT0220 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0, "GSAT0220") # B01
    GSAT0214 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, C_RAAN, 0, "GSAT0214") # C01
    GSAT0226 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, A_RAAN, 0, "GSAT0226", A02_offset) # A02
    #GSAT0221 = Orbits.new_orbit(GalileoRadius + R_Earth, GalileoEcc, GalileoInc, B_RAAN, 0) # B02

    satellites = (GSAT0218, GSAT0220, GSAT0214, GSAT0226)
    
    reception_time = rand(range(0.0, 86400, step=1))
    user_cartesian = [reception_time, 0, 0, R_Earth] # north pole in cartesian coordinates (x,y,z)

    search_duration = 0.2 # in seconds
    step_size = 1e-7

    interval_func = CerealInv.get_interval_func(user_cartesian, GSAT0218)
    user_sat_intervals = CerealInv.get_separation_intervals(interval_func, user_cartesian[1], search_duration, step_size) # return a TimeSeries for the user-satellite separation spacetime intervals
    #search_times = collect(range(user_cartesian[1] - search_duration, user_cartesian[1] + search_duration, step=step_size))

    print("Reception time is: $(user_cartesian[1])")
    plot(user_sat_intervals.times, user_sat_intervals.vals)


end

main()