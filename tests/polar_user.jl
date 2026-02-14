include("../orbits.jl")
include("../kepler.jl")
include("../reltools.jl")

using cereal
using LinearAlgebra
using Plots
import .Orbits
import .RelTools

pyplot()

# structure to hold the results of computing the emission time
# of one satellite with respect to the user location
struct SimulationReport
    reception_time
    emission_time
    orbital_times_set
    orbital_time_step
    separation_intervals
    nullest_interval
    interval_deltas
end

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

# define the interval of a satellite's orbit from which it was physically possible for a signal to have reached the user location by reception time
# the two endpoints of this interval are the furthest time in the past and nearest time in the past from which a signal could have been emitted
function get_orbital_times(reception_time, sat, user)
    c = 3e+8
    # compute all possible travel distances between user location and sat throughout sat's orbit
    distances = []
    for time in range(0.0, 67897.0, step=1)
        distance = R(time, sat, user)
        push!(distances, distance)
    end

    all_travel_durations = distances/c # the corresponding signal travel times for distances[]
    max_duration = maximum(all_travel_durations)
    min_duration = minimum(all_travel_durations)
    durations_range = range(min_duration, max_duration, length=33948) # range of physically plausible signal travel times
    orbital_times_set = [] # where we'll store the portion of the orbital period corresponding to durations_range[]

    # populate orbital_times_range
    for duration in reverse(durations_range)
        push!(orbital_times_set, reception_time - duration)
    end
    
    return orbital_times_set
end

# get all the emission coordinates four-vectors over an interval of orbital timestamps for a given satellite
function get_sat_4vectors(orbital_times_range, sat)
    emission_coords_range = []

    for orbital_time in orbital_times_range
        push!(emission_coords_range, [orbital_time, sat.x(orbital_time), sat.y(orbital_time), sat.z(orbital_time)])
    end

    return emission_coords_range
end

# given set of satellite four-vectors and the user four-vector, compute the corresponding separation four-vectors
function get_sat_user_separations(emission_coords_range, user)
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

# calculate the Δinterval for each time step in a set of intervals
function get_interval_deltas(intervals)
    deltas = []
    prev = nothing
    counter = 1

    for interval in intervals
        if isnothing(prev)
            prev = interval
            counter += 1
            continue
        else
            push!(deltas, (interval - prev))
            prev = interval
            counter += 1
        end
    end

    return deltas
end

# given a set of spacetime intervals, this function returns the interval closest to 0
function find_nullest_interval(interval_pairs)
    nullest_pair = nothing
    for pair in interval_pairs
        if isnothing(nullest_pair)
            nullest_pair = pair
        else
            if norm(pair[2] - 0) < norm(nullest_pair[2] - 0)
                nullest_pair = pair
            end
        end
    end

    return nullest_pair
end

# Run a simulation from generation of reception time to computation of interval sets for one satellite.
# Returns: a fully-populated SatelliteSimulation{} structure.
function run_satellite_simulation(t, sat, user)
    orbital_times = get_orbital_times(t, sat, user)
    orbital_times_range = maximum(orbital_times) - minimum(orbital_times)
    orbital_time_step = orbital_times_range/(length(orbital_times) - 1)

    sat_4vectors = get_sat_4vectors(orbital_times, sat)
    sat_user_separations = get_sat_user_separations(sat_4vectors, [t; user])
    separation_intervals = get_separations_intervals(sat_user_separations)
    interval_deltas = get_interval_deltas(separation_intervals)

    
    interval_pairs = zip(orbital_times, separation_intervals) # an interval pair is a (orbital_time, separation_interval) tuple
    nullest_interval_pair = find_nullest_interval(interval_pairs)

    return SimulationReport(t, nullest_interval_pair[1], orbital_times, orbital_time_step, separation_intervals, nullest_interval_pair, interval_deltas)
    
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

    #=
    # the portions of the satellite's orbits from which a signal have possibly been transmitted to the user
    orbital_times_sets = [
        get_orbital_times(reception_time, GSAT0218, user_cartesian),
        get_orbital_times(reception_time, GSAT0220, user_cartesian),
        get_orbital_times(reception_time, GSAT0214, user_cartesian),
        get_orbital_times(reception_time, GSAT0226, user_cartesian)
    ]

    # the set of spacetime vectors for a satellite along the portion of its orbit obtained in the previous step
    sat_4vector_sets = [
        get_sat_4vectors(orbital_times_sets[1], GSAT0218),
        get_sat_4vectors(orbital_times_sets[2], GSAT0220),
        get_sat_4vectors(orbital_times_sets[3], GSAT0214),
        get_sat_4vectors(orbital_times_sets[4], GSAT0226)
    ]

    # the range of separation four-vectors between the satellite and the user location
    separation_sets = [
        get_sat_user_separations(sat_4vector_sets[1], [reception_time; user_cartesian]),
        get_sat_user_separations(sat_4vector_sets[2], [reception_time; user_cartesian]),
        get_sat_user_separations(sat_4vector_sets[3], [reception_time; user_cartesian]),
        get_sat_user_separations(sat_4vector_sets[4], [reception_time; user_cartesian])
    ]

    # the spacetime intervals of the separation vectors from the previous step
    separation_interval_sets = [
        get_separations_intervals(separation_sets[1]),
        get_separations_intervals(separation_sets[2]),
        get_separations_intervals(separation_sets[3]),
        get_separations_intervals(separation_sets[4])
    ]

    # the per time-step change in spacetime interval for each satellite
    separation_interval_delta_sets = [
        get_interval_deltas(separation_interval_sets[1]),
        get_interval_deltas(separation_interval_sets[2]),
        get_interval_deltas(separation_interval_sets[3]),
        get_interval_deltas(separation_interval_sets[4])
    ]
    =#

    GSAT0218SimulationReport = run_satellite_simulation(reception_time, GSAT0218, user_cartesian)

    

    
    #savefig(plot(orbital_times_ranges[1], separation_interval_sets[1], yaxis="User-Satellite Spacetime Interval", xaxis="Orbital Time (ECI frame)", title="Simulation for reception time=$(reception_time)", legend=false), "img\\North_Pole_Tests\\GSAT0218_test.png")
    


end

main()