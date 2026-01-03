include("orbits.jl")
include("kepler.jl")
include("reltools.jl")

using cereal
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
    
    #
    for i in range(1,10,10)
        ECI_time = rand(0:86400)
        τ1 = GSAT0218.τ(ECI_time)
        τ2 = GSAT0220.τ(ECI_time) + 3600
        τ3 = GSAT0214.τ(ECI_time) + 7400
        τ4 = GSAT0226.τ(ECI_time) + 11400
        #τ5 = GSAT0221.τ(ECI_time)
        
        X1 = [τ1, GSAT0218.x(τ1), GSAT0218.y(τ1), GSAT0218.z(τ1)] 
        X2 = [τ2, GSAT0220.x(τ2), GSAT0220.y(τ2), GSAT0220.z(τ2)]
        X3 = [τ3, GSAT0214.x(τ3), GSAT0214.y(τ3), GSAT0214.z(τ3)]
        X4 = [τ4, GSAT0226.x(τ4), GSAT0226.y(τ4), GSAT0226.z(τ4)]
        X = [X1 X2 X3 X4]
        
        print("ECI time = $(ECI_time)")
        print("\n\n")

        if RelTools.check_constellation(X)
            Xc = locator(X)
            valid_1 = RelTools.check_solution(Xc[1], X)
            valid_2 = RelTools.check_solution(Xc[2], X)

            if Xc[1] == Xc[2]
                if valid_1
                    print("There is only one solution and it is valid\n\n")
                    display(Xc[1])
                    print("\n")
                else
                    print("There is only one 'solution' and it is not valid\n\n")
                    display(Xc[1])
                    print("\n")
                end
            else
                if valid_1 && valid_2
                    print("There are two valid solutions\n\n")
                    display(Xc[1])
                    print("\n")
                    display(Xc[2])
                    print("\n")
                elseif valid_1 || valid_2
                    print("There were two solutions and only this one was valid:\n")
                    if valid_1
                        display(Xc[1])
                        print("\n")
                    else
                        display(Xc[2])
                        print("\n")
                    end
                else
                    print("Neither of these two are valid:\n")
                    display(Xc[1])
                    display(Xc[2])
                    print("\n")
                end
            end
        else
            print("Satellites are not spacelike-separated")
            print("\n\n\n")
        end

    end

    # 10 days in seconds, punctuated by the hour
    #times = range(0.0, 861640.0, step=60)

    #plot(GSAT0218_r.(times))
    #plot([GSAT0218.x.(times), GSAT0220.x.(times), GSAT0214.x.(times)], [GSAT0218.y.(times), GSAT0220.y.(times), GSAT0214.y.(times)], [GSAT0218.z.(times), GSAT0220.z.(times), GSAT0214.z.(times)]) # plotting GSAT0218
    
    
    #plot(times, GSAT0218.τ.(times), xlabel="Time Since Perigee (in seconds)", ylabel="Galileo Satellite's Proper Time")

    #plot(times, Kepler.E.(times, moon.elements.a, moon.elements.ϵ), xlabel="Time Since Perigee (in seconds)", ylabel="Eccentric Anomaly (in degrees)")
    #plot(times, moon.r.(times), xlabel="Time Since Perigee (in seconds)", ylabel="Earth-Moon Distance (in meters)")
    #plot(times, moon.φ.(times), xlabel="Time Since Perigee (in seconds)", ylabel="Moon's Orbital Angle")
    #plot(times, moon.v.(times), xlabel="Time Since Perigee (in seconds)", ylabel="Moon Orbital Velocity (in m/s)")
    
    #savefig("img/moon_radius.png")

end

main()