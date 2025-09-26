include("orbits.jl")
include("kepler.jl")

using .Kepler
using .Orbits
using Plots

#=

Lunar Osculating Orbit Elements:
    a = 384400e+3m
    ϵ = 0.05490
    i = 5.145°
    Ω = 348°
    ω = 67.6°

I'm writing these here for recording purposes, will probably refactor
kepler_constructor() to take an array of arguments soon

=#

function main()
    a = 384400e+3
    ϵ = 0.05490
    i = 5.145
    Ω = 348
    ω = 67.6

    true_anomalies = Array{Float64, 1}()
    true_cosines = Array{Float64, 1}()
    times = Array{Float64, 1}()
    E = 0
    ϕ = nothing # true anomaly, in degrees


    # for each hour over 27 days, calculate the true anomaly
    # this is intended to model the Moon's orbit
    for t in 0.0:3600.0:5184000.0
        append!(times, t)
        E = Kepler.E(E, ϵ, a, t)
        true_cosine = (ϵ - cos(E))/(ϵ*cos(E) - 1)
        append!(true_cosines, true_cosine)
        #ϕ = acos((ϵ - cos(E))/(ϵ*cos(E) - 1))
        #append!(true_anomalies, E)
    end

    r = true_cosine->a*(1-ϵ^2)/(1+ϵ*true_cosine)
    print(length(times))
    print("\n")
    print(length(true_anomalies))
    plot(times, r.(true_cosines), xlabel="Time Since Perigee (in seconds)", ylabel="Earth-Moon Distance (in meters)")
    #savefig("img/moon_radius.png")


    #=print(times)
    print("\n")
    print(ecc_anomalies)=#

end

main()


#=
orbit1 = Orbits.kepler_constructor(a,ϵ,i,Ω,ω)
orbit2 = Orbits.kepler_constructor(a,ϵ,i + 90.0,Ω,ω)


φ = range(0,360)
plot([orbit1.x.(φ), orbit2.x.(φ)], [orbit1.y.(φ), orbit2.y.(φ)], [orbit1.z.(φ), orbit2.z.(φ)], title="Two Orthogonal Orbits", label=["Orbit 1" "Orbit 2"])
savefig("code/img/orthogonal_orbits.png")=#