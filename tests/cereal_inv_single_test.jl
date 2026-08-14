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
    
    A01_emission_vector = CerealInv.get_emission_vector(user_cartesian, GSAT0218)

    print("Emission time is: $(A01_emission_vector[1])\n\n")

    emission_interval = RelTools.get_interval(A01_emission_vector - user_cartesian)

    print("Separation spacetime interval at emission time is: $(emission_interval)\n\n")

    #=interval_func = CerealInv.get_interval_func(user_cartesian, GSAT0226)
    interval_series = CerealInv.get_separation_intervals(interval_func, reception_time, 0.10*GSAT0226.T, 0.001)
    l_bracket, r_bracket = CerealInv.get_root_brackets(interval_series)
    max_time = interval_series.times[indexin(maximum(interval_series.vals), interval_series.vals)[1]]


    positives = sum(x->(x>0), interval_series.vals)
    brackets = find_brackets(interval_series.vals)
    print("Reception time is: $(reception_time)\n\n")
    print("There are $(length(brackets)) crossings and $(positives) positive interval values\n\n")
    print("First crossing is between t=$(interval_series.times[brackets[1][1]]) and t=$(interval_series.times[brackets[1][2]]) | Second crossing is between t=$(interval_series.times[brackets[2][1]]) and t=$(interval_series.times[brackets[2][2]])\n\n")
    #print("The two crossings start at: $(interval_series.times[crossings[1]]) and $(interval_series.times[crossings[2]])\n\n")=#
    
    #plot(interval_series.times, interval_series.vals)

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