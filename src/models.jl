# Export
export slew_time_single_axis
export compute_los_vector
export constraint_agility_single_axis

##########################
# Spacecraft Slew Models #
##########################

"""
Compute the required slew time to maneuver from pointing aligned with the start
vector to the end vector.

Arguments:
- `z_start::Array{<:Real, 1}` Initial pointing vector
- `z_end::Array{<:Real, 1}` Final pointing vector
- `slew_rate::Real` Slew rate [deg/s]

Returns:
- `slew_time::Float64` Time required to slew from start to end axis
"""
function slew_time_single_axis(z_start::Array{<:Real, 1}, z_end::Array{<:Real, 1}; slew_rate=1.0::Real)

    if length(z_start) != 3
        throw(ErrorException("Length of line-of-sight vector z_start must be 3 not $(length(z_start))."))
    end

    if length(z_end) != 3
        throw(ErrorException("Length of line-of-sight vector z_end must be 3 not $(length(z_end))."))
    end

    # Use the dot product to get the angle between two vectors
    proj = dot(z_start, z_end)/(norm(z_start)*norm(z_end))
    
    # Normalize argument if greather than one
    if proj > 1.0
        proj = proj/norm(proj)
    end
    
    # Compute slew angle
    slew_angle = acos(proj)*RAD2DEG

    # Compute slew time
    slew_time = slew_angle/slew_rate

    return slew_time
end 

"""
Compute line of sight vector 
"""
function compute_los_vector(opp::Opportunity, epc::Epoch)
    # Compute starting and ending attitude
    r_eci_location = rECEFtoECI(epc) * opp.location.ecef
    r_eci_sat = state(opp.spacecraft.tle, epc)[1:3]

    # Compute initial look angle
    z_los = r_eci_location - r_eci_sat # Compute look angle 

    # Normalize vector
    z_los = z_los/norm(z_los)

    return z_los
end

"""
Computes whether it is feasible for a spacecraft to slew from the start collect
to the end collect. The transition is based on the spacecraft ability to satisfy
single_axis pointing constraints.

Arguments:
- `start_collect::Opportunity` Initial collect. Spacecraft is assumed pointing here.
- `end_collect::Opportunity` End collect. Spacecraft is assumed pointing here.

Returns:
- `feasible::Bool` `true` if the transition from start to end is feasible. `false` otherwise
"""
function constraint_agility_single_axis(start_collect::Opportunity, end_collect::Opportunity)
    # Check inputs
    if start_collect.spacecraft.id != end_collect.spacecraft.id
        throw(ErrorException("Collect opportunities $(start_collect.id) and $(end_collect.id) for different spacecraft."))
    end

    # Get spacecraft slew-rate
    slew_rate = start_collect.spacecraft.slew_rate

    # Maximum slew time based on worst-case reorigntation from slew-rate
    max_slew_time = 180.0/slew_rate

    # Can't go backwards in time
    if end_collect.t_start < start_collect.t_end
        return false
    end

    # Exit early if time separation is large enough to guarantee feasibility
    if (end_collect.t_start - start_collect.t_end) > max_slew_time
        return true
    end

    # Compute start and end line of sight vectors
    z_start = compute_los_vector(start_collect, start_collect.t_end)
    z_end   = compute_los_vector(end_collect, end_collect.t_start)

    # Compute slew time between orientations
    t_slew = slew_time_single_axis(z_start, z_end, slew_rate=slew_rate)

    # Check if slew time is less than the time available to complete the slew
    if t_slew > (end_collect.t_start - start_collect.t_end)
        return false
    else
        return true
    end
end