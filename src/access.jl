# Exports
export access_geometry, visible
export compute_access
export parallel_compute_access
export create_lookup_location_opportunity
export create_lookup_location
export create_lookup_opportunity

###################
# Access Geometry #
###################

"""
Compute the view geometry from an observer to a specific request. 

Arguments:
- `sat_ecef:Array{<:Real, 1}` Earth-fixed satellite position
- `request::Request` Request object

Returns:
- `look_angle::Float64` Look angle from the satellite to the request center. Equivalent to off-nadiar angle. [deg]
- `range::Float64` Range from the observing satellite to the location. [m]
"""
function access_geometry(sat_ecef::Array{<:Real, 1}, request::Request)
    # Satellite state
    r_sat = sat_ecef[1:3]
    
    # Geodetic sub-satellte point
    sat_geod 	 = sECEFtoGEOD(r_sat)
    sub_sat_geod = [sat_geod[1], sat_geod[2], 0.0]
    sub_sat_ecef = sGEODtoECEF(sub_sat_geod)
    
    # Look angle
    nadir_dir  = (sub_sat_ecef - r_sat)/norm(sub_sat_ecef - r_sat)
    location_dir = (request.ecef - r_sat)/norm(request.ecef - r_sat)
    look_angle = acos(dot(nadir_dir, location_dir))*180.0/pi

    # Distance from sub-satellite point to location along direct line of sight
    range = norm(r_sat - request.ecef) # range to location

    return look_angle, range
end

"""
Compute the view geometry from an observer to a specific station. 

Arguments:
- `sat_ecef:Array{<:Real, 1}` Earth-fixed satellite position
- `station::GroundStation` Station object

Returns:
- `elevation::Float64` Elevation of satellite with respect to station [deg]
"""
function access_geometry(sat_ecef::Array{<:Real, 1}, station::GroundStation)
    # Satellite state
    r_sat     = sat_ecef[1:3]
    r_station = station.ecef

    azimuth, elevation, range = sSEZtoAZEL(sECEFtoSEZ(r_station, r_sat), use_degrees=true)

    return elevation, range
end


"""
Computes whether an request is visible from a given state and return boolean true/false
result.

Arguments:
- `sat_ecef::Array{<:Real, 1}` Satellite position in Earth fixed frame
- `request::Request` Request being viewed

Returns:
- `visible::Bool` Indication of whether request is visible or not
"""
function visible(sat_ecef::Array{<:Real, 1}, request::Request)
    # Compute access geometry
    look_angle, range = access_geometry(sat_ecef, request)

    r = norm(sat_ecef[1:3])
    # theta_max = request.look_angle_max*pi/180.0
    # max_range = -r*cos(theta_max) - sqrt(r^2*cos(theta_max)^2 - (r^2-R_EARTH^2)^2)
    
    # To be fixed, no idea why this is complex:
    max_range = R_EARTH 
    # max_range = 1.5*sqrt(r^2 - R_EARTH^2)

    # @debug "max range: $max_range"
    if request.look_angle_min <= look_angle &&
       look_angle <= request.look_angle_max &&
       range <= max_range

        return true
    else
        return false
    end
end

"""
Computes whether an station is visible from a given state and return boolean true/false
result.

Arguments:
- `sat_ecef::Array{<:Real, 1}` Satellite position in Earth fixed frame
- `station::GroundStation` Station being communicatted with

Returns:
- `visible::Bool` Indication of whether station is visible or not
"""
function visible(sat_ecef::Array{<:Real, 1}, station::GroundStation)
    elevation, range = access_geometry(sat_ecef, station)

    r = norm(sat_ecef[1:3])

    if station.elevation_min <= elevation
        # dist = norm(sat_ecef[1:3] - station.ecef[1:3])
        return true
    else
        return false
    end
end


###################
# Location Access #
###################

"""
Find the time when a location transitions between being accessible/not-accessible.
"""
function find_access_boundary(tle::TLE, location::Location, epc0::Epoch, step::Real; tol::Real=1.0e-1)

    # println("Step size: $step")

    epc = deepcopy(epc0)

    if abs(step) < tol
        # If time step is below desired tolerance do nothing
        return epc0
    else
        x_ecef   = sECItoECEF(epc, state(tle, epc))
        initial_visibility = visible(x_ecef, location)

        # While visibility is the same as it was initially step in time and 
        # recompute visibility
        while visible(x_ecef, location) == initial_visibility
            epc    = epc + step
            x_ecef = sECItoECEF(epc, state(tle, epc))
        end

        # Reverse search direction and half step size
        next_step = -sign(step)*max(abs(step)/2, tol/2)

        return find_access_boundary(tle, location, epc, next_step, tol=tol)
    end
end

function compute_access(tle::TLE, locations::Array{<:Location, 1}, 
            epc_min::Union{Epoch, Nothing}, epc_max::Union{Epoch, Nothing}; 
            macro_step::Real=60.0, tol::Real=0.01, id_offset::Integer=0)
    if epc_min == nothing
        epc_min = tle.epoch
    end

    if epc_max == nothing
        epc_min += 86400.0
    end

    # Compute orbital period
    T = orbit_period(sCARTtoOSC(state(tle, epc_min), use_degrees=true)[1])

    opportunities = Array{Opportunity, 1}(undef, 0)

    next_step = macro_step

    for loc in locations
        # Perform macro adaptive stepsize search
        epc = epc_min
        while epc < epc_max

            # println("Current step: $epc")

            x_ecef = sECItoECEF(epc, state(tle, epc))

            if visible(x_ecef, loc)
                # println("Found instant of visibility - $(visible(x_ecef, loc)) - $epc")

                # Search for AOS (before initial guess epoch)
                window_open = find_access_boundary(tle, loc, epc, -macro_step, tol=tol)
                # println("Found window_open boundary: $window_open")

                # Search for LOS (after initial guess epoch)
                window_close = find_access_boundary(tle, loc, epc, macro_step, tol=tol)
                # println("Found window_close boundary: $window_close")

                teststep = (window_close - window_open)*1.1

                # If zero-doppler collection is required updated pass-times and geometry profile
                collect_duration = 0
                if typeof(loc) == Request && loc.require_zero_doppler == true
                    epc_mid = window_open + (window_close - window_open)/2.0
                    window_open  = epc_mid - loc.collect_duration/2.0
                    window_close = epc_mid + loc.collect_duration/2.0
                end

                # Create Oppourtunity
                if (window_close - window_open) > 0.0
                    if typeof(loc) == Request
                        collect = Collect(window_open, 
                            window_close,
                            orbit=tle,
                            location=loc)

                        push!(opportunities, collect)
                    elseif typeof(loc) == GroundStation
                        contact = Contact(window_open, 
                            window_close,
                            orbit=tle,
                            location=loc)

                        push!(opportunities, contact)
                    else
                        throw(ErrorException("Unknown type of location: $loc"))
                    end
                end

                # Step a most of an orbit because there (likely) won't be another
                # Acces until at least 1 orbit later
                # epc += T*(13/16)
                epc += T*(1/2)
                # epc += T/10
                # epc += teststep

            else
                # Take a macro step because there's a good chance there won't be an
                # opportunity until the next one
                epc += macro_step
            end
        end
    end

    # Sort opportunities in ascending order
    sort!(opportunities, by = x -> x.t_start)

    # Apply IDs for collects in ascending time order
    for (idx, opp) in enumerate(opportunities)
        # Override computed opportunity id
        opp.id = idx + id_offset
    end

    return opportunities
end

function parallel_compute_access(tle::TLE, locations::Array{<:Location, 1},
            epc_min::Union{Epoch, Nothing}, epc_max::Union{Epoch, Nothing}; 
            macro_step::Real=60.0, tol::Real=0.01, id_offset::Integer=0)

    # Create anonymous function to map array inputs to 
    fn(a, b, c, d, e) = x -> compute_access(a, x, b, c, macro_step=d, tol=e)

    # Create work assignments
    nw = nworkers()        # Number of workers
    lw = length(locations) # Length of work

    # @debug "Have $lw items of work and $nw workers")

    wa = floor(Int, lw/nw) # Average work per worker
    wr = lw - nw*wa        # Remaining worker

    # Work Array
    assignments = Array{<:Location, 1}[]

    # @debug "Worker 1 convering $(1):$(wa+wr)")
    push!(assignments, locations[1:(wa+wr)])
    for i in 2:nw
        # @debug "Worker $i convering $(1+wr+(i-1)*wa):$(i*wa+wr)")
        push!(assignments, locations[(1+wr+(i-1)*wa):(i*wa+wr)])
    end

    results = pmap(fn(tle, epc_min, epc_max, macro_step, tol), assignments)
    
    # Aggregate all results
    opportunities = vcat(results...)

    # Sort opportunities in ascending order
    sort!(opportunities, by = x -> x.t_start)

    # Aggregate results and resolve id numbering
    for (idx, opp) in enumerate(opportunities)
        # Override computed opportunity id
        opp.id = idx + id_offset
    end

    return opportunities
end

"""
Create dictionary lookup of the possible opportunities that exist for each request.

Arguments:
- `opportunities::Array{<:Opportunity, 1}` Array of opportunities for all requests

Returns:
- `location_access::Dict{<:Location, Array{<:Opportunity, 1}}` Lookup table which returns the array of opportunities for each request
"""
function create_lookup_location_opportunity(locations::Array{<:Location, 1}, opportunities::Array{<:Opportunity, 1})
    # Get ID Types
    loc_id_type = typeof(locations[1].id)
    opp_id_type = typeof(opportunities[1].id)

    # Initialize dictionary storing collects for each request
    location_oppoturnities = Dict{loc_id_type, Array{opp_id_type, 1}}()


    # Create collect array for each unique request
    for loc in locations
        location_oppoturnities[loc.id] = Integer[]
    end

    # Populate lookup with opportunities
    for opp in opportunities
        push!(location_oppoturnities[opp.location.id], opp.id)
    end

    return location_oppoturnities
end

"""
Create Lookup to get location based on location ID
"""
function create_lookup_location(locations::Array{<:Location})
    loc_id_type = typeof(locations[1].id)

    lookup = Dict{loc_id_type, Location}()

    for location in locations
        lookup[location.id] = location
    end

    return lookup
end

"""
Create Lookup to get opportunity based on opportunity ID
"""
function create_lookup_opportunity(opportunities::Array{<:Opportunity})
    loc_id_type = typeof(opportunities[1].id)

    lookup = Dict{loc_id_type, Opportunity}()

    for opportunity in opportunities
        lookup[opportunity.id] = opportunity
    end

    return lookup
end

