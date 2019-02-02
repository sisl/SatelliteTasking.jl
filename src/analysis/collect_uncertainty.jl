__precompile__(true)
module CollectUncertainty

# Julia Imports
using Statistics
using SatelliteDynamics.Time: Epoch

# Package Imports
using SatelliteTasking.DataStructures: Orbit, Image, Collect
using SatelliteTasking.Collection: find_all_collects, collect_diff, collect_stats, find_matching_collect
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
    miss_img = Int64[] # Store number of missing images per window bin

    for i in 1:num_bins
        # Hourly bins of windows
        epc_window_start = epc_min + 3600*(i-1)
        epc_window_end   = epc_min + 3600*i

        @debug "Computing statistings for window $(epc_window_start) - $(epc_window_end)"
        
        # TODO: Refactor the use of collect stats here to be efficient. This duplicates a bit of work
        win_mean, win_sdev, win_miss = collect_stats(true_collects, all_collects, epc_min=epc_window_start, epc_max=epc_window_end)

        push!(mean_diff, [win_mean...])
        push!(sdev_diff, [win_sdev...])
        push!(miss_img, win_miss)
    end

    # Arrange into matrix
    mean_diff = hcat(mean_diff...)
    sdev_diff = hcat(sdev_diff...)

    return true_collects, perturbed_collects, mean_diff, sdev_diff, miss_img
end

export find_missing_collect
"""
Find the collects that are different between two sets of collects.

Arguments:
- `collect_list_a` First list of collect
- `collect_list_b` Second list of collects

Returns
- `missing_from_a` List of collects present in list b but _not_ present in a
- `missing_from_b` List of collects present in list a but _not_ present in b
"""
function find_missing_collect(collect_list_a::Array{Collect, 1}, collect_list_b::Array{Collect, 1})
    # Initialize missing lists
    missing_from_a = Collect[]
    missing_from_b = Collect[]

    # First check for collects from b missing in a
    for col_b in collect_list_b
        if find_matching_collect(collect_list_a, col_b) == nothing
            push!(missing_from_a, col_b)
        end
    end

    # Second check for collects in a missing from b
    for col_a in collect_list_a
        if find_matching_collect(collect_list_b, col_a) == nothing
            push!(missing_from_b, col_a)
        end
    end

    return missing_from_a, missing_from_b
end

end # End CollectUncertainty module
