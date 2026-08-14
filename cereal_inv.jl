module CerealInv

include("reltools.jl")

using DoubleFloats
using Roots
import .RelTools

export generate_emission_constellation

struct TimeSeries
    times
    vals
end

function generate_emission_constellation(user, sats)
    emission_vectors = []
    
    for sat in sats
        emission_vector = get_emission_vector(user, sat)
        push!(emission_vectors, emission_vector)
    end

    return emission_vectors
end

function get_emission_vector(user, sat)
    interval_func = get_interval_func(user, sat)
    search_duration = Double64(0.2) # in seconds
    step_size = Double64(1e-7) # in seconds
    user_sat_intervals = get_separation_intervals(interval_func, user[1], search_duration, step_size) # return a TimeSeries for the user-satellite separation spacetime intervals

    #max_index = indexin(maximum(user_sat_intervals.vals), interval_series.vals)[1]
    #max_time = interval_series.times[max_index]
    l_index, r_index = get_root_brackets(user_sat_intervals.vals)[1]
    l_bracket, r_bracket = user_sat_intervals.times[l_index], user_sat_intervals.times[r_index]
    emission_time = find_emission_time(interval_func, l_bracket, r_bracket) # find the null interval using some root-finding function

    return [emission_time, sat.x(emission_time), sat.y(emission_time), sat.z(emission_time)]
end

function get_interval_func(user, sat)
    return t->(RelTools.get_interval([t, sat.x(t), sat.y(t), sat.z(t)] - user))
end

# evaluate user-satellite intervals over a time range centered on the reception time and with a specified step step_size
function get_separation_intervals(interval_func, reception_time, search_duration, step_size)
    times = collect(range(reception_time - search_duration, reception_time, step=step_size))
    intervals = []

    for time in times
        push!(intervals, interval_func(time))
    end

    return TimeSeries(times, intervals)
end

# search through a set of user-satellite separation intervals and store
# tuples containing the indices of pairs of intervals with opposite signs
function get_root_brackets(intervals)
    i = 0
    brackets = []
    for interval in intervals
        i += 1
        if i != length(intervals)
            if intervals[i+1] * interval < Double64(0)
                push!(brackets, (i, i+1))
            end
        end
    end

    return brackets
end

# apply an algorithm from Roots.jl using
# the tuple of indices representing a bracket
function find_emission_time(interval_func, l_bracket, r_bracket)
    return find_zero(interval_func, (l_bracket, r_bracket), Roots.A42())
end

end