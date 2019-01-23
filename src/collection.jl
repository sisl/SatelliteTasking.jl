__precompile__(true)
module Collection

# Julia imports
using LinearAlgebra
using SatelliteDynamics.Constants: R_EARTH
using SatelliteDynamics.Coordinates: sECEFtoGEOD, sGEODtoECEF

# Package imports
using SatelliteTasking.DataStructures: Orbit, Image, Collect

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

    @debug "max range: $max_range"

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
    run = []
    groups = [run]
    expected_idx = nothing

    for idx in indices
        if expected_idx == nothing || idx == expected_idx
            push!(groups[end], idx)
        else
            run = []
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
- `collects::Array{Collect, 1}`
"""
function groups_extract_collects(orbit::Orbit, image::Image, index_groups)
    # collects = Array{Collect, 1}(undef, length(index_groups))
    collects = []
    
    for (i, group) in enumerate(index_groups)
        if length(group) > 0
            sow = orbit.epc[group[1]]
            eow = orbit.epc[group[end]]
            # collects[i] = Collect(sow, eow, id=-1, orbit_id=orbit.id, image_id=image.id)
            push!(collects, Collect(sow, eow, id=-1, orbit_id=orbit.id, image_id=image.id))
        end
    end
    
    return collects
end

export find_collects
"""
Find all collects a given image is visible for an orbit.

Arguments:
- `orbit::Orbit` Orbit of observing satellite
- `image::Image` Image under observation

Returns:
- `collects::Array{Collect, 1}` Array of collection colortunities given the orbit
"""
function find_collects(orbit::Orbit, image::Image)
    visibility = Array{Bool, 1}(undef, length(orbit.t))

    for i in 1:length(orbit.t)
        visibility[i] = image_visible(orbit.ecef[:, i], image)
    end

    visible_indices = group_indices(findall(visibility))
    collects   = groups_extract_collects(orbit, image, visible_indices)

    return collects
end

export find_all_collects
"""
Find all collection colortunities for a set of images.

Arguments:
- `orbit::Orbit` Orbit of observing satellite
- `images::Array{Image, 1}` Array of images 
- `sort::Bool` Sort collects in ascending order by start of collect window

Returns:
- `collects::Array{Collect, 1}` Array of collection colortunities for all images
"""
function find_all_collects(orbit::Orbit, images::Array{Image, 1}; sort=true::Bool)
    collects = []
    for img in images
        for col in find_collects(orbit, img)
            push!(collects, col)
        end
    end

    # Sort collects in start-of-window order
    if sort
        sort!(collects, by = x -> x.sow)
    end

    return collects
end


"""
Internal helper function. Return the midtime for a collect.

Arguments:
- `collect::Collect` Collect to compute midtime for

Returns:
- `epc::Epoch` Midtime collect window
"""
function collect_midtime(collect::Collect)
    return collect.sow + (collect.eow - collect.sow)/2.0
end

"""
Finds the closest matching collect out of a list 

Arguments:
- `collect_list` List of collections to extract closest match from
- `col::Collect` Collect to match from list

Returns:
- `collect::Collect` The closest matching in `collect_list` to the input collect.
"""
function find_matching_collect(collect_list::Array{<:Any, 1}, col::Collect)
    col_midtime = collect_midtime(col)
    for o in collect_list
        if o.sow < col_midtime && o.eow > col_midtime && o.image_id == col.image_id
           return o
        end
    end

    # If no matching collect found return nothing
    return nothing
end

export collect_diff
"""
Compute the difference in collection times between two different sets of collections.
Finds matching sets of collects between the two collect lists. 

Arguments:
- `collect_list_a::Array{Collect, 1}` First set of collections
- `collect_list_b::Array{Collect, 1}` Second set of collections

Returns:
- `collect_diffs::Array{Array{Float64, 1}, 1}` Array of differences between collect windows. Returns difference in start time, end time, and duration for each matched collect between the two lists.
"""
function collect_diff(collect_list_a::Array{<:Any, 1}, collect_list_b::Array{<:Any, 1})
    col_diffs = []
    for col in collect_list_b
        matching_a = find_matching_collect(collect_list_a, col)
        if matching_a != nothing
            sow_diff = col.sow - matching_a.sow
            eow_diff = col.eow - matching_a.eow
            dur_diff = col.duration - matching_a.duration
            push!(col_diffs, [sow_diff, eow_diff, dur_diff])
        end
    end

    return col_diffs
end

export collect_pair
"""
Gets the subset of Collects common to two sets of Collects.

Arguments:
- `collect_list_a::Array{Collect, 1}` First set of collections
- `collect_list_b::Array{Collect, 1}` Second set of collections

Returns:
- `collect_diffs::Array{Collect, 1}` Array of Collects common to the two lists
"""
function collect_pair(collect_list_a::Array{<:Any, 1}, collect_list_b::Array{<:Any, 1})
    col_pairs = []
    for col in collect_list_b
        matching_a = find_matching_collect(collect_list_a, col)
        if matching_a != nothing
            push!(col_pairs, (matching_a, col))
        end
    end

    return col_pairs
end

export collect_stats
"""
Compute the statistics on the differences of the collect windows of collect list
B with respect to collect list A.

Arguments:
- `collect_list_a::Array{Collect, 1}` First set of collections
- `collect_list_b::Array{Collect, 1}` Second set of collections

Returns:
- `collect_stats::Tuple{Array{Float64, 1}, Array{Float64, 1}}`: Mean and standard deviation of the differences in start time, end time, and collection window duration.
"""
function collect_stats(collect_list_a::Array{<:Any, 1}, collect_list_b::Array{<:Any, 1}; epc_min=nothing, epc_max=nothing)
    collects = []

    # Extract Opportunities from each list
    for cols in collect_list_b
        c = cols
        if epc_min != nothing
            c = filter(x -> x.sow > epc_min, c)
        end

        if epc_max != nothing
            c = filter(x -> x.sow < epc_max, c)
        end
        push!(collects, c)
    end

    collects = vcat(collects...)
    
    # Compute the difference between all collects and the list of "true" collects
    col_diffs = collect_diff(collect_list_a, collects)

    # Concatenate the matrix into a single matrix to easily compute statistic
    col_errors = hcat(col_diffs...)

    # Compute mean and standard deviation
    col_mean = mean(col_errors, dims=2)
    col_sdev = std(col_errors, dims=2)

    return col_mean, col_sdev
end

end # End module