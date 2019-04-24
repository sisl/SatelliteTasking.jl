__precompile__(true)
module Collection

# Julia imports
using LinearAlgebra
using Statistics

# SatelliteDynamics imports
using SatelliteDynamics.Constants: R_EARTH
using SatelliteDynamics.Coordinates: sECEFtoGEOD, sGEODtoECEF, sECEFtoENZ, sENZtoAZEL

# Package imports
using SatelliteTasking.DataStructures: Orbit, Image, Location, GroundStation, Opportunity

using SatelliteDynamics


export image_view_geometry
"""
Compute the view geometry from an observer to a specific image. 

Arguments:
- `sat_ecef:Array{<:Real, 1}` Earth-fixed satellite position
- `image::Image` Image object

Returns:
- `look_angle::Float64` Look angle from the satellite to the image center. Equivalent to off-nadiar angle. [deg]
- `range::Float64` Range from the observing satellite to the location. [m]
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
    location_dir = (image.ecef - r_sat)/norm(image.ecef - r_sat)
    look_angle = acos(dot(nadir_dir, location_dir))*180.0/pi

    # Distance from sub-satellite point to location along direct line of sight
    range = norm(r_sat - image.ecef) # range to location

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
    max_range = 1.5*sqrt(r^2 - R_EARTH^2)

    # @debug "max range: $max_range"
    if image.look_angle_min <= look_angle &&
       look_angle <= image.look_angle_max &&
       eow_range <= max_range

        return true
    else
        return false
    end
end

export station_view_geometry
"""
Compute the view geometry from an observer to a specific station. 

Arguments:
- `sat_ecef:Array{<:Real, 1}` Earth-fixed satellite position
- `station::GroundStation` Station object

Returns:
- `elevation::Float64` Elevation of satellite with respect to station [deg]
"""
function station_view_geometry(sat_ecef::Array{<:Real, 1}, station::GroundStation)
    # Satellite state
    r_sat     = sat_ecef[1:3]
    r_station = station.ecef

    azimuth, elevation, range = sSEZtoAZEL(sECEFtoSEZ(r_station, r_sat), use_degrees=true)

    return elevation, range
end

export station_visible
"""
Computes whether an station is visible from a given state and return boolean true/false
result.

Arguments:
- `sat_ecef::Array{<:Real, 1}` Satellite position in Earth fixed frame
- `station::GroundStation` Station being communicatted with

Returns:
- `visible::Bool` Indication of whether station is visible or not
"""
function station_visible(sat_ecef::Array{<:Real, 1}, station::GroundStation, time)
    elevation, range = station_view_geometry(sat_ecef, station)

    r = norm(sat_ecef[1:3])

    if station.elevation_min <= elevation
        dist = norm(sat_ecef[1:3] - station.ecef[1:3])
        # println("$time - $elevation - $range - $(dist/1e3) - $(sECEFtoGEOD(sat_ecef, use_degrees=true))")

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
function groups_extract_opportunities(orbit::Orbit, location::Location, index_groups)
    # opportunities = Array{Opportunity, 1}(undef, length(index_groups))
    opportunities = Opportunity[]
    
    for (i, group) in enumerate(index_groups)
        if length(group) > 0
            sow = orbit.epc[group[1]]
            eow = orbit.epc[group[end]]
            if typeof(location) == Image
                push!(opportunities, Opportunity(sow, eow, orbit=orbit, location=location, collect_duration=location.collect_duration))
            else
                push!(opportunities, Opportunity(sow, eow, orbit=orbit, location=location))
            end
        end
    end
    
    return opportunities
end

export find_opportunities
"""
Find all opportunities a given location is visible for an orbit.

Arguments:
- `orbit::Orbit` Orbit of observing satellite
- `location::Location` Image under observation

Returns:
- `opportunities::Array{Opportunity, 1}` Array of collection opportunities given the orbit
"""
function find_opportunities(orbit::Orbit, location::Location)
    visibility = Array{Bool, 1}(undef, length(orbit.t))

    for i in 1:length(orbit.t)
        if typeof(location) == Image
            visibility[i] = image_visible(orbit.ecef[:, i], location)
        elseif typeof(location) == GroundStation
            visibility[i] = station_visible(orbit.ecef[:, i], location, orbit.epc[i])
        end
    end

    visible_indices = group_indices(findall(visibility))
    opportunities   = groups_extract_opportunities(orbit, location, visible_indices)

    return opportunities
end

export find_all_opportunities
"""
Find all collection opportunities for a set of locations.

Arguments:
- `orbit::Orbit` Orbit of observing satellite
- `locations::Array{Location, 1}` Array of locations 
- `sort::Bool` Sort opportunities in ascending order by start of window

Returns:
- `opportunities::Array{Opportunity, 1}` Array of collection opportunities for all locations
"""
function find_all_opportunities(orbit::Orbit, locations::Array{<:Location, 1}; sort::Bool=true, zero_doppler::Bool=false)
    opportunities = Opportunity[]
    for tar in locations
        for opp in find_opportunities(orbit, tar)
            push!(opportunities, opp)
        end
    end

    # If zero-doppler constrained edit all opportunities start and and end times
    # to border the mid time.
    if zero_doppler == true
        for opp in opportunities
            opp.sow = opp.mid - opp.location.collect_duration/2.0
            opp.eow = opp.mid + opp.location.collect_duration/2.0
            opp.duration = opp.eow - opp.sow
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
    opportunity_candidates = filter(x -> x.location == col.location, opportunity_list)

    # Screen each candidate and return if below match threshold
    for c in opportunity_candidates
        if abs(c.mid - col.mid) < match_max_seconds
            return c
        end
    end

    @debug "Missing opportunity: $col"

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

    for opp_b in opportunity_list_b
        matching_a = find_matching_opportunity(opportunity_list_a, opp_b)
        if matching_a != nothing
            sow_diff = opp_b.sow - matching_a.sow
            eow_diff = opp_b.eow - matching_a.eow
            dur_diff = opp_b.duration - matching_a.duration

            push!(opp_diffs, [sow_diff, eow_diff, dur_diff])
        else
            @debug "Unable to find matching opportunity for $opp_b"
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
function opportunity_stats(opportunity_list_a::Array{Opportunity, 1}, opportunity_list_b::Array{Opportunity, 1})    
    # Compute the difference between all opportunities and the list of "true" opportunities
    opp_diffs, opp_miss = opportunity_diff(opportunity_list_a, opportunity_list_b)

    # Concatenate the matrix into a single matrix to easily compute statistic
    opp_errors = hcat(opp_diffs...)

    # Compute mean and standard deviation
    opp_mean, opp_sdev = zeros(Float64, 3), zeros(Float64, 3)
    if length(opp_errors) > 0
        opp_mean = mean(opp_errors, dims=2)
        opp_sdev = std(opp_errors, dims=2)
    end

    return opp_mean, opp_sdev, opp_miss
end

export split_opportunities
"""
Compute discrete collection opportunities for 

Arguments:
- `opportunity_list::Array{Opportunity,1}` List of opportunities to compute collections for
- `max_opps::Integer` Maximum number of collects to divide each opportunity into

Returns:
- `collects::Array{Opportunity,1}` Array of collection opportunities
"""
function split_opportunities(opportunity_list::Array{Opportunity,1}, max_opps=0::Integer)
    opps = Opportunity[]

    # Compute Opportunity for each opportunity
    for opportunity in opportunity_list
        next_opp_start = deepcopy(opportunity.sow)
        num_opps       = 0

        if max_opps == 0
            max_opps = floor((opportunity.eow - opportunity.sow)/opportunity.collect_duration)
        end

        # Stupider algorithm for dividing up collect windows
        while num_opps < max_opps && (next_opp_start+opportunity.collect_duration) < opportunity.eow
            push!(opps, Opportunity(next_opp_start, 
                                    next_opp_start+opportunity.collect_duration,
                                    orbit=opportunity.orbit,
                                    location=opportunity.location))

            next_opp_start += opportunity.collect_duration
            num_opps       += 1
        end
    end
    
    # Sort opps so they are in ascending order
    sort!(opps, by = x -> x.sow)
    
    return opps
end

export group_image_opportunities
"""
Create dictionary lookup of the possible opportunities that exist for each image.

Arguments:
- `opportunities::Array{Collect, 1}` Array of opportunities for all images

Returns:
- `image_opportunities::Dict{Image, Array{Collect, 1}}` Lookup table which returns the array of opportunities for each image
"""
function group_image_opportunities(opportunities::Array{Opportunity, 1})
    # Initialize dictionary storing collects for each image
    image_opportunities = Dict{Image, Array{Tuple{Int64, Opportunity}, 1}}()

    # Create collect array for each unique image
    for opp in opportunities
        image_opportunities[opp.location] = Tuple{Int64, Opportunity}[]
    end

    # Populate lookup with opportunities
    for (i, opp) in enumerate(opportunities)
        push!(image_opportunities[opp.location], (i, opp))
    end

    return image_opportunities
end

end # End module