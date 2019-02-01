__precompile__(true)
module CollectUncertainty

# Julia Imports
using Statistics
using SatelliteDynamics.Time: Epoch

# Package Imports
using SatelliteTasking.DataStructures: Orbit, Image, Collect
using SatelliteTasking.Collection: find_all_collects, collect_diff, collect_stats
using SatelliteTasking.Simulation: simulate_orbits

export compute_perturbed_collects
"""
Computes the true collection opportunities, as well as the collection opporrtunities
for the perturbed orbits.

All statistics are computed for the following collect prorperties (sow, eow, diff)

Arguments:
- `true_orbit::Orbit` True orbit object used to compute the collections considered "the truth"
- `perturbed_orbits::Array{Orbit,1}` List of perturbed orbits
- `images::Array{Images,1}` Set of images for which to compute and analyze collect times
- `epc_min::Epoch` _Optional_ Start of time window to compute statistics for. If none is provided the entire orbit arc will be used.
- `epc_max::Epoch` _Optional_ End of time window to compute statistics for. If none is provided the entire orbit arc will be used.
- `epc_step::Epoch` _Optional_ Time step used to bin statistics analysis. Expects units of seconds

Returns:
- `true_collects::Array{Collect,1}` List of collect collection opportunities
- `perturbed_collects::Arrray{Array{Collect,1}}` List of collects computed for each perrturbed orbit 
- `mean_diff` Mean of differernces in collection properties between perturbed and true orbits. 
- `sdev_diff` Standard deviation of differernces in collection properties between perturbed and true orbits.
"""
function compute_perturbed_collects(true_orbit::Orbit, perturbed_orbits::Array{Orbit,1}, images::Array{Image,1}; epc_min=nothing, epc_max=nothing, epc_step=nothing)
    
    # Validate inputs 
    if epc_min != nothing && !(true_orbit.epc[1] <= epc_min <= true_orbit.epc[end])
        throw(ArgumentError("epc_min outside of orbit timespan"))
    end

    if epc_max != nothing && !(true_orbit.epc[1] <= epc_max <= true_orbit.epc[end])
        throw(ArgumentError("epc_max outside of orbit timespan"))
    end

    # Extract number of orbits from input
    num_orbits = length(perturbed_orbits)

    # Compute true collects
    true_collects = find_all_collects(true_orbit, images, sort=true)

    # Compute collections
    perturbed_collects = Array{Collect, 1}[]
    collect_diffs      = Array{Float64, 1}[]
    for i in 1:num_orbits
        push!(perturbed_collects, find_all_collects(perturbed_orbits[i], images, sort=true))
    end

    # Aggregate collections form all perturbed orbits
    all_collects = sort!(vcat(perturbed_collects...), by = x -> x.sow)

    # Compute difference in windows between initial and final
    collect_diffs    = collect_diff(true_collects, all_collects)
    collect_diff_mat = hcat(collect_diffs...)
    
    # Compute statistics on collect differences
    mean_diff = Array{Float64, 1}[]
    sdev_diff = Array{Float64, 1}[]

    # Initialize Statistics windows is none is provided
    if epc_min == nothing
        epc_min = true_orbit.epc[1]
    end

    if epc_max == nothing
        epc_max = true_orbit.epc[end]
    end

    if epc_step == nothing
        epc_step = epc_max - epc_min
    end

    # Validate step as last part
    if epc_step != nothing && !(mod(epc_max-epc_min, epc_step) == 0.0)
        throw(ArgumentError("epc_step does not evenly divide time span"))
    end

    # Compute statistics for each leveling bin
    num_bins = Int((epc_max - epc_min)/epc_step)

    for i in 1:num_bins
        # Hourly bins of windows
        epc_window_start = epc_min + 3600*(i-1)
        epc_window_end   = epc_min + 3600*i

        @debug "Computing statistings for window $(epc_window_start) - $(epc_window_end)"
        
        # TODO: Refactor the use of collect stats here to be efficient. This duplicates a bit of work
        m, s = collect_stats(true_collects, [all_collects], epc_min=epc_window_start, epc_max=epc_window_end)

        push!(mean_diff, [m...])
        push!(sdev_diff, [s...])
    end

    # Arrange into matrix
    mean_diff = hcat(mean_diff...)
    sdev_diff = hcat(sdev_diff...)

    return true_collects, perturbed_collects, mean_diff, sdev_diff 
end

end # End CollectUncertainty module
