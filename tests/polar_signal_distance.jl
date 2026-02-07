include("../orbits.jl")
include("../kepler.jl")
include("../reltools.jl")

using cereal
using LinearAlgebra
using Plots
import .Orbits
import .RelTools

pyplot()

# returns the Euclidean distance between the satellite and user at time t
function R(t, sat, user)
    X = [sat.x(t), sat.y(t), sat.z(t)]
    S = X - user
    return sqrt(S[1]^2 + S[2]^2 + S[3]^2)
end

function get_emission_points(emission_times, sat)
    emission_points = []
    for emission_time in emission_times
        push!(emission_points, [emission_time, sat.x(emission_time), sat.y(emission_time), sat.z(emission_time)])
    end

    return emission_points
end

# select a physically realistic satellite position that is temporally prior to reception time at north pole
function get_orbital_times_range(reception_time, sat, user)
    c = 3e+8
    # compute all possible travel distances between user location and GSAT0218
    distances = []
    for time in range(0.0, 67897.0, step=1)
        distance = R(time, sat, user)
        push!(distances, distance)
    end

    all_travel_durations = distances/c # the corresponding signal travel times for distances[]
    max_duration = maximum(all_travel_durations)
    min_duration = minimum(all_travel_durations)
    durations_range = range(min_duration, max_duration, length=33948) # range of physically plausible signal travel times
    orbital_times_range = [] # where we'll store the portion of the orbital period corresponding to durations_range[]

    # populate orbital_times_range
    for duration in reverse(durations_range)
        push!(orbital_times_range, reception_time - duration)
    end
    
    return orbital_times_range
end

# get all the emission coordinates four-vectors over an interval of orbital timestamps for a given satellite
function get_emission_coords_range(orbital_times_range, sat)
    emission_coords_range = []

    for orbital_time in orbital_times_range
        push!(emission_coords_range, [orbital_time, sat.x(orbital_time), sat.y(orbital_time), sat.z(orbital_time)])
    end

    return emission_coords_range
end

# given set of emission coordinates and a user location, compute the corresponding separation vectors
function get_emission_user_separations(emission_coords_range, user)
    separations = []

    for emission_coords in emission_coords_range
        push!(separations, emission_coords - user)
    end

    return separations
end

# given set of separation vectors, compute their respective spactime intervals
function get_separations_intervals(separations)
    c = 3e+8
    intervals = []
    
    for separation in separations
        push!(intervals, -(c*separation[1])^2 + separation[2]^2 + separation[3]^2 + separation[4]^2)
    end

    return intervals
end

# given a set of spacetime intervals, this function returns the interval closest to 0
function find_nullest_interval(intervals)
    nullest = intervals[1]
    for interval in intervals
        if norm(interval - 0) < norm(nullest - 0)
            nullest = interval
        end
    end

    return nullest
end

function next_emission_coord(emission_coord, sat)
    new_emission_time = emission_coord[1] + 1

    return [new_emission_time, sat.x(new_emission_time), sat.y(new_emission_time), sat.z(new_emission_time)]
end

function find_true_emission_coord(emission_coord, user, sat)
    S = [emission_coord[1] - user[1], emission_coord[2] - user[2], emission_coord[3] - user[3], emission_coord[4] - user[4]]

    if RelTools.is_lightlike(S)
        return emission_coord
    else
        return find_true_emission_coord(next_emission_coord(emission_coord, sat), user, sat)
    end
end

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
    
    user_spherical = [R_Earth, 0, 0] # north pole in polar coordinates (r,θ,φ)
    user_cartesian = [0, 0, R_Earth] # north pole in cartesian coordinates (x,y,z)
    reception_time = rand(range(0, 86400, step=1)) # reception time will be randomly sampled from the first 24 hours of orbit

    orbital_times_ranges = [
        get_orbital_times_range(reception_time, GSAT0218, user_cartesian),
        get_orbital_times_range(reception_time, GSAT0220, user_cartesian),
        get_orbital_times_range(reception_time, GSAT0214, user_cartesian),
        get_orbital_times_range(reception_time, GSAT0226, user_cartesian)
    ]

    emission_coords_ranges = [
        get_emission_coords_range(orbital_times_ranges[1], GSAT0218),
        get_emission_coords_range(orbital_times_ranges[2], GSAT0220),
        get_emission_coords_range(orbital_times_ranges[3], GSAT0214),
        get_emission_coords_range(orbital_times_ranges[4], GSAT0226)
    ]

    separations = [
        get_emission_user_separations(emission_coords_ranges[1], [reception_time; user_cartesian]),
        get_emission_user_separations(emission_coords_ranges[2], [reception_time; user_cartesian]),
        get_emission_user_separations(emission_coords_ranges[3], [reception_time; user_cartesian]),
        get_emission_user_separations(emission_coords_ranges[4], [reception_time; user_cartesian])
    ]

    separation_interval_sets = [
        get_separations_intervals(separations[1]),
        get_separations_intervals(separations[2]),
        get_separations_intervals(separations[3]),
        get_separations_intervals(separations[4])
    ]


    p = plot([orbital_times_ranges[1], orbital_times_ranges[2], orbital_times_ranges[3], orbital_times_ranges[4]], [separation_interval_sets[1], separation_interval_sets[2], separation_interval_sets[3], separation_interval_sets[4]], label=["GSAT0218" "GSAT0220" "GSAT0214" "GSAT0226"], layout=(4,1))
    display(p)
    sleep(30)
    
    #print("Reception time: $(reception_time)\n\n")
    


end

main()