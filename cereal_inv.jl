module CerealInv

using DoubleFloats
import FiniteDifferences

export find_sat

# structure to hold the results of computing the emission time
# of one satellite with respect to the user location
struct SimulationResult
    reception_time
    emission_time
    orbital_times_set
    orbital_time_step
    separation_intervals
    nullest_interval
    interval_deltas
    interval_derivs
    position
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

function zoom(nullest_point, sat, user)
    c = 3e+8
    m = 1.9e+16
    delta_interval = nullest_point[2]
    time_step = delta_interval/m
    new_time = nothing
    
    if nullest_point[2] < 0
        new_time = nullest_point[1] + time_step
    elseif nullest_point[2] > 0
        new_time = nullest_point[1] - time_step
    end

    new_sat_X = [new_time, sat.x(new_time), sat.y(new_time), sat.z(new_time)]
    S = new_sat_X - user
    return -(c*S[1])^2 + S[2]^2 + S[3]^2 + S[4]^2
end






# define the interval of a satellite's orbit from which it was physically possible for a signal to have reached the user location by reception time
# the two endpoints of this interval are the furthest time in the past and nearest time in the past from which a signal could have been emitted
function get_orbital_times(reception_time, sat, user, full_orbit=false)
    c = 3e+8

    # compute all possible travel distances between user location and sat throughout sat's orbit
    distances = []
    for time in range(0.0, 67897.0, step=1)
        distance = R(time, sat, user)
        push!(distances, distance)
    end

    all_travel_durations = distances/c # the corresponding signal travel times for distances[]
    orbital_times_set = []

    if full_orbit
        orbital_times_set = collect(range(0.0, 67897.0, step=1)) # all times at which the satellite position is computed
    else
        max_duration = maximum(all_travel_durations)
        min_duration = minimum(all_travel_durations)
        durations_range = range(min_duration, max_duration, length=33948) # range of physically plausible signal travel times
        #orbital_times_set = [] # where we'll store the portion of the orbital period corresponding to durations_range[]

        # populate orbital_times_range
        for duration in reverse(durations_range)
            push!(orbital_times_set, reception_time - duration)
        end
    end

    return orbital_times_set
end

# get all the emission coordinates four-vectors over an interval of orbital timestamps for a given satellite
function get_sat_4vectors(orbital_times, sat)
    sat_4vectors = []

    for orbital_time in orbital_times
        push!(sat_4vectors, [orbital_time, sat.x(orbital_time), sat.y(orbital_time), sat.z(orbital_time)])
    end

    return sat_4vectors
end

# given set of satellite four-vectors and the user four-vector, compute the corresponding separation four-vectors
function get_sat_user_separations(sat_4vectors, user)
    separations = []

    for sat_4vector in sat_4vectors
        push!(separations, sat_4vector - user)
    end

    return separations
end

# return spacetime interval of any 4-vector
function get_spacetime_interval(vector)
    c = 3e+8

    return -(c*vector[1])^2 + vector[2]^2 + vector[3]^2 + vector[4]^2
end

# given set of separation vectors, compute their respective spactime intervals
function get_separations_intervals(separations)
    intervals = []
    
    for separation in separations
        push!(intervals, get_spacetime_interval(separation))
    end

    return intervals
end

# the full function that get_sat_interval_func
# calls behind the scenes
function get_separation_interval(t, sat, user)
    sat_4vector = [t, sat.x(t), sat.y(t), sat.z(t)]
    
    separation = sat_4vector - [t; user]

    interval = get_spacetime_interval(separation)

    return interval
end

# takes full context (satellite AND user) and returns
# anonymous one-dimensional function of time for calculation
# of spacetime interval between specified sat & user
function get_sat_interval_func(sat, user)
    return t->get_separation_interval(t, sat, user)
end

function get_interval_derivs(interval_func, times, time_step)
    # Finite difference functions
    m_central = FiniteDifferences.central_fdm(5, 1)  
    m_forward = FiniteDifferences.forward_fdm(5, 1)
    m_backward = FiniteDifferences.backward_fdm(5, 1)
    
    # Create container for derivatives
    derivs = similar(times)

    for i in eachindex(times)
        if i == 1
            # Forward difference at the start boundary
            derivs[i] = m_forward(interval_func, times[i])
        elseif i == length(times)
            # Backward difference at the end boundary
            derivs[i] = m_backward(interval_func, times[i])
        else
            # Central difference for the interior
            derivs[i] = m_central(interval_func, times[i])
        end
    end

    return derivs
end

function get_significant_digits(num)
    sig = significand(num)
    sig_string = string(sig)
    sig_figs = length(sig_string) - 1

    return sig_figs
end

function report_roundings(nums)
    rounding_count = 0

    for num in nums
        sigs = get_significant_digits(num)
        if sigs == 16
            rounding_count += 1
        end
    end

    return rounding_count

end

# Run a simulation from generation of reception time to computation of interval sets for one satellite.
# Returns: a fully-populated SatelliteSimulation{} structure.
function find_sat(t, sat, user)
    orbital_times = get_orbital_times(t, sat, user, true)
    orbital_times_range = maximum(orbital_times) - minimum(orbital_times)
    #orbital_time_step = orbital_times_range/(length(orbital_times) - 1)
    orbital_time_step = 1

    sat_4vectors = get_sat_4vectors(orbital_times, sat)
    sat_user_separations = get_sat_user_separations(sat_4vectors, [t; user])
    
    separation_intervals = get_separations_intervals(sat_user_separations)
    intervals_roundings = report_roundings(separation_intervals)
    intervals_roundings_percentage = (intervals_roundings / length(separation_intervals)) * 100
    print("Interval Separations has $(intervals_roundings) cases of 16-digit significands, which is $(intervals_roundings_percentage)% of all separation intervals.\n\n")
    
    #interval_deltas = get_interval_deltas(separation_intervals)
    interval_derivs = get_interval_derivs(get_sat_interval_func(sat, user), orbital_times, orbital_time_step)
    derivs_roundings = report_roundings(interval_derivs)
    derivs_roundings_percentage = (derivs_roundings / length(interval_derivs)) * 100
    print("Interval Derivative has $(derivs_roundings) cases of 16-digit significands, which is $(derivs_roundings_percentage)% of all interval derivatives.\n\n")

    
    interval_pairs = zip(orbital_times, separation_intervals) # an interval pair is a (orbital_time, separation_interval) tuple
    #nullest_intervaltime = find_nullest_interval(interval_pairs)

    pos = [sat.x(t), sat.y(t), sat.z(t)]
    
    return SimulationResult(t, nothing, orbital_times, orbital_time_step, separation_intervals, nothing, nothing, interval_derivs, pos)
    
end

end