__precompile__(true)
module CollectUncertainty

# Julia Imports
using Statistics
using SatelliteDynamics.Time: Epoch

# Package Imports
using SatelliteTasking.DataStructures: Orbit, Image, Opportunity, Collect
using SatelliteTasking.Collection: find_all_opportunities, opportunity_diff, opportunity_stats, find_matching_opportunity, compute_collects_by_number
using SatelliteTasking.Simulation: simulate_orbits

export compute_perturbed_opportunities
"""
Computes the true collection opportunities, as well as the collection opportunities
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
- `true_opportunities::Array{Opportunity,1}` List of collect collection opportunities
- `perturbed_opportunities::Arrray{Array{Opportunity,1}}` List of opportunities computed for each perrturbed orbit 
- `mean_diff` Mean of differernces in opporunity properties between perturbed and true orbits. 
- `sdev_diff` Standard deviation of differernces in opportunities properties between perturbed and true orbits.
"""
function compute_perturbed_opportunities(true_orbit::Orbit, perturbed_orbits::Array{Orbit,1}, images::Array{Image,1}; epc_min=nothing, epc_max=nothing, epc_step=nothing)
    
    # Validate inputs 
    if epc_min != nothing && !(true_orbit.epc[1] <= epc_min <= true_orbit.epc[end])
        throw(ArgumentError("epc_min outside of orbit timespan"))
    end

    if epc_max != nothing && !(true_orbit.epc[1] <= epc_max <= true_orbit.epc[end])
        throw(ArgumentError("epc_max outside of orbit timespan"))
    end

    # Extract number of orbits from input
    num_orbits = length(perturbed_orbits)

    # Compute true opportunities
    true_opportunities = find_all_opportunities(true_orbit, images, sort=true)

    # Compute collections
    perturbed_opportunities = Array{Opportunity, 1}[]
    opportunity_diffs       = Array{Float64, 1}[]
    for i in 1:num_orbits
        @debug "Computing opportunities for perturbed orbit: $i"
        push!(perturbed_opportunities, find_all_opportunities(perturbed_orbits[i], images, sort=true))
    end

    # Aggregate collections form all perturbed orbits
    all_opportunities = sort!(vcat(perturbed_opportunities...), by = x -> x.sow)
    
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

        # Extract Opportunities from each list
        window_opportunities = filter(x -> epc_window_start <= x.sow < epc_window_end, all_opportunities)

        @debug "Found $(length(window_opportunities)) opportunities in window"
        
        # TODO: Refactor the use of collect stats here to be efficient. This duplicates a bit of work
        win_mean, win_sdev, win_miss = opportunity_stats(true_opportunities, window_opportunities)

        push!(mean_diff, [win_mean...])
        push!(sdev_diff, [win_sdev...])
        push!(miss_img, win_miss)
    end

    # Arrange into matrix
    mean_diff = hcat(mean_diff...)
    sdev_diff = hcat(sdev_diff...)

    return true_opportunities, perturbed_opportunities, mean_diff, sdev_diff, miss_img
end

# export compute_perturbed_collects
# """
# Computes true and pertrubed collect opportunities to be used in planning 
# """
# function compute_perturbed_collects(true_opportunities::Array{Opportunity, 1}, perturbed_opportunities::Array{Opportunity, 1}; max_collects=0::integer)
# end

export find_missing_opportunity
"""
Find the opportunities that are different between two sets of opportunities.

Arguments:
- `opportunity_list_a` First list of collect
- `opportunity_list_b` Second list of opportunities

Returns
- `missing_from_a` List of opportunities present in list b but _not_ present in a
- `missing_from_b` List of opportunities present in list a but _not_ present in b
"""
function find_missing_opportunity(opportunity_list_a::Array{Opportunity, 1}, opportunity_list_b::Array{Opportunity, 1})
    # Initialize missing lists
    missing_from_a = Opportunity[]
    missing_from_b = Opportunity[]

    @debug opportunity_list_a
    @debug opportunity_list_b

    # First check for opportunities from b missing in a
    for opp_b in opportunity_list_b
        if find_matching_opportunity(opportunity_list_a, opp_b) == nothing
            push!(missing_from_a, opp_b)
        end
    end

    # Second check for opportunities in a missing from b
    for opp_a in opportunity_list_a
        if find_matching_opportunity(opportunity_list_b, opp_a) == nothing
            push!(missing_from_b, opp_a)
        end
    end

    return missing_from_a, missing_from_b
end

end # End CollectUncertainty module
