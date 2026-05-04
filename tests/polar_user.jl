include("../orbits.jl")
include("../kepler.jl")
include("../reltools.jl")
include("../cereal_inv.jl")

using cereal
using LinearAlgebra
using Plots
using .CerealInv
import .Orbits
import .RelTools

#pyplot()

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
    
    user_spherical = [R_Earth, 0, 0] # north pole in polar coordinates (r,θ,φ)
    user_cartesian = [0, 0, R_Earth] # north pole in cartesian coordinates (x,y,z)
    

    

    test_number = 1
    test_dir_root = "C:\\FengLab\\code\\kepler_simulator\\results\\North_Pole_Tests\\Series21"

    for i in range(1,6)
        #reception_time = rand(range(0, 86400, step=1)) # reception time will be randomly sampled from the first 24 hours of orbit
        reception_time = rand(range(0.0, 67897.0, step=1))
        test_dir_path = test_dir_root * "\\Constellation$(test_number)"
        mkdir(test_dir_path)

        for sat in satellites
            result = find_sat(reception_time, sat, user_cartesian)

            sat_test_path = test_dir_path * "\\$(sat.name)"
            mkdir(sat_test_path)

            savefig(plot(result.orbital_times_set, result.separation_intervals, title="$(sat.name) Test $(test_number): Interval vs. Time", legend=false, yaxis="Sat-User Separation Interval (m^2)", xaxis="Time in Seconds"), sat_test_path * "\\interval_graph.png")
            #savefig(plot(result.orbital_times_set, result.interval_deltas, title="$(sat.name) Test $(test_number): Interval Deltas vs. Time", legend=false, yaxis="Separation Interval Deltas (m^2)", xaxis="Time in Seconds"), sat_test_path * "\\delta_graph.png")
            savefig(plot(result.orbital_times_set, result.interval_derivs, title="$(sat.name) Test $(test_number): Interval Derivative vs. Time", legend=false, yaxis="Separation Interval Derivative (m^2/s)", xaxis="Time in Seconds"), sat_test_path * "\\deriv_graph.png")
            
            test_data = "Reception Time: $(result.reception_time)\nEstimated Emission Time: $(result.emission_time)\nSimulation Timestep Size: $(result.orbital_time_step)\nNullest Interval: result.nullest_interval\nSatellite Position: $(result.position)\n"
            write(sat_test_path * "\\data.txt", test_data)
        end
        
        test_number += 1
    end

end

main()