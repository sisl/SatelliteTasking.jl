__precompile__(true)
module Collection

# Julia imports
using LinearAlgebra
using Statistics

# SatelliteDynamics imports
using SatelliteDynamics.Constants: R_EARTH
using SatelliteDynamics.Coordinates: sECEFtoGEOD, sGEODtoECEF

# Package imports
using SatelliteTasking.DataStructures: Orbit, Image, Opportunity, Collect

export image_view_geometry
"""
Compute the view geometry from an observer to a specific image. 

Arguments:
- `sat_ecef:Array{<:Real, 1}` Earth-fixed satellite position
- `image::Image` Image object

Returns:
- `look_angle::Float64` Look angle from the satellite to the image center. Equivalent to off-nadiar angle. [deg]
- `range::Float64` Range from the observing satellite to the target. [m]
"""
function image_view_geometry(sat_ecef::Array{<:Real, 1}, image::Image)
    # Satellite state
    r_sat = sat_ecef[1:3]
    
    # Geodetic sub-satellte point
    sat_geod 	 = sECEFtoGEOD(r_sat)
    sub_sat_geod = [sat_geod[1], sat_geod[2], 0.0]
    sub_sat_ecef = sGEODtoECEF(sub_sat_geod)
    
    # Look angle
    nadir_dir  = (sub_sat_ecef - r_sat)/norm(sub_sat_ecef - r_sat)
    target_dir = (image.ecef - r_sat)/norm(image.ecef - r_sat)
    look_angle = acos(dot(nadir_dir, target_dir))*180.0/pi

    # Distance from sub-satellite point to target along direct line of sight
    range = norm(r_sat - image.ecef) # range to target

    return look_angle, range
end


export image_visible
"""
Computes whether an image is visible from a given state and return boolean true/false
result.

Arguments:
- `sat_ecef::Array{<:Real, 1}` Satellite position in Earth fixed frame
- `image::Image` Image being viewed

Returns:
- `visible::Bool` Indication of whether image is visible or not
"""
function image_visible(sat_ecef::Array{<:Real, 1}, image::Image)
    look_angle, eow_range = image_view_geometry(sat_ecef, image)

    r         = norm(sat_ecef[1:3])
    # theta_max = image.look_angle_max*pi/180.0
    # max_range = -r*cos(theta_max) - sqrt(r^2*cos(theta_max)^2 - (r^2-R_EARTH^2)^2)
    
    # To be fixed, no idea why this is complex:
    max_range = sqrt(r^2 - R_EARTH^2)

    # @debug "max range: $max_range"

    if image.look_angle_min <= look_angle &&
       look_angle <= image.look_angle_max &&
       eow_range <= max_range
        return true
    else
        return false
    end
end

"""
Internal helper function to group an array of linear increasing indices 
indicating periods of visibility.

Arguments:
- `indices::Array{<:Int, 1}` List of indices

Returns:
- `groups::Array{Array{Int32, 1}, 1}` Array of arrays where each sub-array contains a continuous set of indicies where the image is visible.
"""
function group_indices(indices::Array{<:Int, 1})
    run = Int64[]
    groups = Array{Int64,1}[run]
    expected_idx = nothing

    for idx in indices

        if expected_idx == nothing || idx == expected_idx
            # Add current index to run
            push!(groups[end], idx)
        else
            # Start new run if necessary
            run = Int64[idx]
            push!(groups, run)
        end

        expected_idx = idx + 1
    end

    return groups
end

"""
Internal helper function to create 
Arguments:
- `orbit::Orbit` Orbit 

Returns:
- `opportunities::Array{Opportunity, 1}`
"""
function groups_extract_opportunities(orbit::Orbit, image::Image, index_groups)
    # opportunities = Array{Opportunity, 1}(undef, length(index_groups))
    opportunities = Opportunity[]
    
    for (i, group) in enumerate(index_groups)
        if length(group) > 0
            sow = orbit.epc[group[1]]
            eow = orbit.epc[group[end]]
            push!(opportunities, Opportunity(sow, eow, orbit=orbit, 
                                    image=image, dwell_time=image.dwell_time))
        end
    end
    
    return opportunities
end

export find_opportunities
"""
Find all opportunities a given image is visible for an orbit.

Arguments:
- `orbit::Orbit` Orbit of observing satellite
- `image::Image` Image under observation

Returns:
- `opportunities::Array{Opportunity, 1}` Array of collection opportunities given the orbit
"""
function find_opportunities(orbit::Orbit, image::Image)
    visibility = Array{Bool, 1}(undef, length(orbit.t))

    for i in 1:length(orbit.t)
        visibility[i] = image_visible(orbit.ecef[:, i], image)
    end

    visible_indices = group_indices(findall(visibility))
    opportunities   = groups_extract_opportunities(orbit, image, visible_indices)

    return opportunities
end

export find_all_opportunities
"""
Find all collection opportunities for a set of images.

Arguments:
- `orbit::Orbit` Orbit of observing satellite
- `images::Array{Image, 1}` Array of images 
- `sort::Bool` Sort opportunities in ascending order by start of window

Returns:
- `opportunities::Array{Opportunity, 1}` Array of collection opportunities for all images
"""
function find_all_opportunities(orbit::Orbit, images::Array{Image, 1}; sort=true::Bool)
    opportunities = Opportunity[]
    for img in images
        for col in find_opportunities(orbit, img)
            push!(opportunities, col)
        end
    end

    # Sort opportunities in start-of-window order
    if sort
        sort!(opportunities, by = x -> x.sow)
    end

    return opportunities
end

"""
Finds the closest matching opportunity out of a list 

Arguments:
- `opportunity_list` List of opportunities to extract closest match from
- `col::Opportunity` Opportunity to match from list
- `match_max_seconds` _Optional_ Maximum difference in seconds to be considered a valid match. Default: 10 min

Returns:
- `opportunity::Opportunity` The closest matching in `opportunity_list` to the input opportunity.
"""
function find_matching_opportunity(opportunity_list::Array{Opportunity, 1}, col::Opportunity; match_max_seconds=600::Real)
    opportunity_candidates = filter(x -> x.image == col.image, opportunity_list)

    # Screen each candidate and return if below match threshold
    for c in opportunity_candidates
        if abs(c.mid - col.mid) < match_max_seconds
            return c
        end
    end

    @debug "Missing opportunity: $col"

    # Older method for matching opportunities
    # for o in opportunity_list
    #     if o.sow < col.mid && o.eow > col.mid && o.image == col.image
    #        return o
    #     end
    # end

    # If no matching opportunity found return nothing
    return nothing
end

export opportunity_diff
"""
Compute the difference in opportunity times between two different sets of opportunities.
Finds matching sets of opportunities between the two opportunity lists. 

Arguments:
- `opportunity_list_a::Array{Opportunity, 1}` First set of opportunities
- `opportunity_list_b::Array{Opportunity, 1}` Second set of opportunities

Returns:
- `opportunity_diffs::Array{Array{Float64, 1}, 1}` Array of differences between opportunity windows. Returns difference in start time, end time, and duration for each matched opportunity between the two lists.
- `opp_miss::Int64` Opportunitys present in b missing from a
"""
function opportunity_diff(opportunity_list_a::Array{Opportunity, 1}, opportunity_list_b::Array{Opportunity, 1})
    opp_diffs = Array{Float64, 1}[]
    opp_miss  = 0
    for opp_a in opportunity_list_a
        matching_b = find_matching_opportunity(opportunity_list_b, opp_a)
        if matching_b != nothing
            sow_diff = matching_b.sow - opp_a.sow
            eow_diff = matching_b.eow - opp_a.eow
            dur_diff = matching_b.duration - opp_a.duration
            push!(opp_diffs, [sow_diff, eow_diff, dur_diff])
        else
            @debug "Unable to find matching opportunity for $opp_a"
            opp_miss += 1
        end
    end

    return opp_diffs, opp_miss
end

export opportunity_stats
"""
Compute the statistics on the differences of the opportunity windows of opportunity list
B with respect to opportunity list A.

Arguments:
- `opportunity_list_a::Array{Opportunity, 1}` First set of opportunities
- `opportunity_list_b::Array{Opportunity, 1}` Second set of opportunities

Returns:
- `opportunity_stats::Tuple{Array{Float64, 1}, Array{Float64, 1}}` Mean and standard deviation of the differences in start time, end time, and collection window duration.
- `opportunity_miss::Int` Number of opportunities missing from opportunity_list_b that are epected to be present in a
"""
function opportunity_stats(opportunity_list_a::Array{Opportunity, 1}, opportunity_list_b::Array{Opportunity, 1}; epc_min=nothing, epc_max=nothing)
    opportunities = copy(opportunity_list_a) # Filter on true opportunities in case perturbed move around

    # Extract Opportunities from each list
    if epc_min != nothing
        opportunities = filter(x -> x.sow > epc_min, opportunities)
    end

    if epc_max != nothing
        opportunities = filter(x -> x.sow < epc_max, opportunities)
    end

    @debug "Found $(length(opportunities)) opportunities in window"
    
    # Compute the difference between all opportunities and the list of "true" opportunities
    opp_diffs, opp_miss = opportunity_diff(opportunities, opportunity_list_b)

    # Concatenate the matrix into a single matrix to easily compute statistic
    opp_errors = hcat(opp_diffs...)

    # Compute mean and standard deviation
    opp_mean, opp_sdev = nothing, nothing
    if length(opp_errors) > 0
        opp_mean = mean(opp_errors, dims=2)
        opp_sdev = std(opp_errors, dims=2)
    end

    return opp_mean, opp_sdev, opp_miss
end

export compute_collects_by_number
"""
Compute discrete collection opportunities for 

Arguments:
- `opportunity_list::Array{Opportunity,1}` List of opportunities to compute collections for
- `max_collects::Integer` Maximum number of collects to divide each opportunity into

Returns:
- `collects::Array{Collects,1}` Array of collection opportunities
"""
function compute_collects_by_number(opportunity_list::Array{Opportunity,1}, max_collects=0::Integer)
    collects = Collect[]

    # Compute Collects for each opportunity
    for opportunity in opportunity_list
        # If maximum number of collects per opportunity is not set, set as the 
        # maximum possible of non-overlapping collects over the opportunity window
        if max_collects == 0
            max_collects = convert(typeof(max_collects), floor(opportunity.duration/opportunity.dwell_time))
        end

        # Time in opportunity window not taken up by a collect
        empty_span = opportunity.duration - max_collects*opportunity.dwell_time

        # Spacing between collects
        collect_spacing = empty_span/(max_collects + 1)

        # Start time of next collect window
        next_collect_start = deepcopy(opportunity.sow)

        # Create collects for opportunities
        for i in 1:max_collects
            next_collect_start += collect_spacing
            push!(collects, Collect(next_collect_start, 
                                    next_collect_start+opportunity.dwell_time,
                                    orbit=opportunity.orbit,
                                    image=opportunity.image,
                                    opportunity=opportunity))
        end
    end
    
    # Sort collects so they are in ascending order
    sort!(collects, by = x -> x.sow)
    
    return collects
end

end # End module