include("orbits.jl")

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
a = 384400e+3
ϵ = 0.05490
i = 5.145
Ω = 348
ω = 67.6


orbit1 = Orbits.kepler_constructor(a,ϵ,i,Ω,ω)
orbit2 = Orbits.kepler_constructor(a,ϵ,i + 90.0,Ω,ω)


φ = range(0,360)
plot([orbit1.x.(φ), orbit2.x.(φ)], [orbit1.y.(φ), orbit2.y.(φ)], [orbit1.z.(φ), orbit2.z.(φ)], title="Two Orthogonal Orbits", label=["Orbit 1" "Orbit 2"])

savefig("code/img/orthogonal_orbits.png")