# Package exports
export view_geometry

"""
Compute the view geometyr from an observer to a specific image.

Arguments:
- `request::Request` Request object
- `scef::Array{<:Rreal, 1}` Earth-fixed satellite position

Returns:
- `look_angle::Real` Off-nadir look angle from satellite to request center
- `range::Real` Range from satellite to location
"""
function view_geometry(request::Request, secef::Array{<:Real, 1})
    # Satellite state
    r_sat = secef[1:3]

    # Request ECEF state
    recef = request.ecef
    
    # Geodetic sub-satellte point
    sat_geod 	 = sECEFtoGEOD(r_sat)
    sub_sat_geod = [sat_geod[1], sat_geod[2], 0.0]
    sub_sat_ecef = sGEODtoECEF(sub_sat_geod)
    
    # Look angle
    nadir_dir  = (sub_sat_ecef - r_sat)/norm(sub_sat_ecef - r_sat)
    location_dir = (request.ecef - r_sat)/norm(request.ecef - r_sat)
    look_angle = acos(dot(nadir_dir, location_dir))*180.0/pi

    # Distance from sub-satellite point to location along direct line of sight
    range = norm(r_sat - image.ecef) # range to location

    return look_angle, range
end

function view_geometry(station::Station, secef::Array{<:Real, 1})
    # Satellite state
    r_sat     = secef[1:3]
    r_station = station.ecef

    azimuth, elevation, range = sSEZtoAZEL(sECEFtoSEZ(r_station, r_sat), use_degrees=true)

    return elevation, range
end


function visible(secef::Array{<:Real, 1}, request::Request)
    look_angle, eow_range = view_geometry(request, secef)

    r = norm(secef[1:3])
    # theta_max = image.look_angle_max*pi/180.0
    # max_range = -r*cos(theta_max) - sqrt(r^2*cos(theta_max)^2 - (r^2-R_EARTH^2)^2)
    
    # To be fixed, no idea why this is complex:
    max_range = 1.5*sqrt(r^2 - R_EARTH^2)

    # @debug "max range: $max_range"
    if request.look_angle_min <= look_angle &&
       look_angle <= request.look_angle_max &&
       eow_range <= max_range

        return true
    else
        return false
    end
end

function visible(station::Station, secef::Array{<:Real, 1})
    elevation, range = view_geometry(station, secef)

    r = norm(secef[1:3])

    if station.elevation_min <= elevation
        dist = norm(secef[1:3] - station.ecef[1:3])

        return true
    else
        return false
    end
end

function visibility_find_boundary(tle::TLE, location::Geolocation, epc0::Epoch, step::Real; tol::Real=1.0e-1)
    """Compute the boundary of when the geolocation 
    """

    # println("Step size: $step")

    epc = deepcopy(epc0)

    if abs(step) < tol
        return epc0
    else
        x_ecef  = sECItoECEF(epc, state(tle, epc))
        ivis = visible(x_ecef, location)

        # println("Initial visibility: $ivis - $epc")

        while visible(x_ecef, location) == ivis
            epc    = epc + step
            x_ecef = sECItoECEF(epc, state(tle, epc))
            # println("Step visibility: $(visible(x_ecef, location)) - $epc")
        end

        # println("Final visibility: $(visible(x_ecef, location)) - $epc")

        next_step = -sign(step)*max(abs(step)/2, tol/2)

        return visibility_find_boundary(tle, location, epc, next_step, tol=tol)
    end
end

function find_all_opportunities(mdp::MDPProblem, location::Array{<:Geolocation, 1}; 
            macro_step::Real=60.0, tol::Real=0.01, zero_doppler::Bool=true)
    
    tle = mdp.tle
    epc_min = mdp.t_start
    epc_max = mdp.t_end
    
    if epc_min == nothing
        epc_min = tle.epoch
    end

    if epc_max == nothing
        epc_min += 86400.0
    end

    # Compute orbital period
    T = orbit_period(sCARTtoOSC(state(tle, epc_min), use_degrees=true)[1])

    opportunities = Opportunity[]

    next_step = macro_step

    for loc in location
        # Perform macro adaptive stepsize search
        epc = epc_min
        while epc < epc_max

            # println("Current step: $epc")

            x_ecef = sECItoECEF(epc, state(tle, epc))

            if visible(x_ecef, loc)
                # println("Found instant of visibility - $(visible(x_ecef, loc)) - $epc")

                # Search for AOS (before initial guess epoch)
                window_open = visibility_find_boundary(tle, loc, epc, -macro_step, tol=tol)
                # println("Found window_open boundary: $window_open")

                # Search for LOS (after initial guess epoch)
                window_close = visibility_find_boundary(tle, loc, epc, macro_step, tol=tol)
                # println("Found window_close boundary: $window_close")

                # If zero-doppler collection is required updated pass-times and geometry profile
                collect_duration = 0
                if (typeof(loc) == Request) && ((window_close - window_open) > 0.0)
                    epc_mid = window_open + (window_close - window_open)/2.0
                    window_open  = get_rel_time(mdp, epc_mid) - loc.duration/2.0
                    window_close = get_rel_time(mdp, epc_mid) + loc.duration/2.0
                    
                    col = Collect(
                        id = length(opportunities) + 1,
                        geolocation_id = loc.id,
                        duration = loc.duration,
                        window_open=window_open,
                        window_close=window_close,
                        reward=loc.reward
                    )
                    push!(opportunities, col)               
                elseif (typeof(loc) == Request) && ((window_close - window_open) > 0.0)
                    epc_mid = window_open + (window_close - window_open)/2.0
                    window_open  = get_rel_time(mdp, epc_mid) - loc.duration/2.0
                    window_close = get_rel_time(mdp, epc_mid) + loc.duration/2.0
                    
                    con = Contact(
                        id = length(opportunities) + 1,
                        geolocation_id = loc.id,
                        duration = loc.duration,
                        window_open=window_open,
                        window_close=window_close
                    )
                    push!(opportunities, con)  
                end

                # Step a half orbit 
                epc += T/2.0

            else
                # Take a macro step because there's a good chance there won't be an
                # opportunity until the next one
                epc += macro_step
            end
        end
    end

    # Sort opportunities in ascending order
    sort!(opportunities, by = x -> x.sow)

    return opportunities
end