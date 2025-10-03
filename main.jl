include("orbits.jl")
include("kepler.jl")

import .Orbits
using Plots

#=

Lunar Osculating Orbit Elements:
    a = 384400e+3m
    ϵ = 0.05490
    i = 5.145°
    Ω = 348°
    ω = 67.6°

=#

function main()
    
    # arguments are: a, ϵ, i, Ω, ω
    moon = Orbits.new_orbit(384400e+3, 0.05490, 5.145, 348, 67.6)


    # for each hour over 60 days, calculate the true anomaly
    # this is intended to model the Moon's orbit
    times = range(0.0, 5184000.0, step=3600)

    #plot(times, Kepler.E.(times, moon.elements.a, moon.elements.ϵ), xlabel="Time Since Perigee (in seconds)", ylabel="Eccentric Anomaly (in degrees)")
    #plot(times, moon.r.(times), xlabel="Time Since Perigee (in seconds)", ylabel="Earth-Moon Distance (in meters)")
    plot(times, moon.φ.(times), xlabel="Time Since Perigee (in seconds)", ylabel="Moon's Orbital Angle")
    #plot(times, moon.v.(times), xlabel="Time Since Perigee (in seconds)", ylabel="Moon Orbital Velocity (in m/s)")
    
    #savefig("img/moon_radius.png")


    #=print(times)
    print("\n")
    print(ecc_anomalies)=#

end

main()