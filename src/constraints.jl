__precompile__(true)
module Constraints

using LinearAlgebra
using SatelliteDynamics: EPoch

using SatelliteTasking.DataStructures: Collect

##########################
# Spacecraft Slew Models #
##########################

"""
Compute 
"""
function slew_time_single_axis(z_start::Array{<:Real, 1}, z_end::Array{<:Real, 1})

    # Use the dot product to get the angle between two vectors
    proj = dot(z_start, z_end)/(norm(z_start)*norm(z_end))
    
    # Normalize argument if greather than one
    if proj > 1.0
        proj = proj/norm(proj)
    end
    
    # Compute slew angle
    slew_angle = acos(proj)*180.0/pi

    return slew_angle
end 

#####################
# Constraint Models #
#####################

# """
# Compute line of sight vector 
# """
# function compute_los_vector(col::Collect, epc::Epoch):
    
#     # Compute starting and ending attitude
#     r_eci_target  = rot_ecef_to_eci(epc) @ col.image.ecef
#     r_eci_sat     = _np.array(col.sat.tle.state(epc))[0:3]

#     # Compute initial look angle
#     z_los = r_eci_target - r_eci_sat

#     # Normalize vector
#     z_los = z_los/_np.linalg.norm(z_los)

#     return z_los

# export constraint_single_axis_slew
# function constraint_single_axis_slew(start_collect::Collect, end_collect::Collect)
# end

# def constraint_single_axis_slew(opp_s, opp_e):
#     '''Computes whether it is feasible to transition from the start opportunity
#     to the end opportunity, given that the transition must satisfy single-axis
#     pointing constraints.

#     Returns True if the transition is feasible, returns false otherwise.
#     '''

#     # Can't go backwards in time
#     if opp_e.aos < opp_s.los:
#         return False

#     # Exit early if time separation is large enough to guarantee feasibility
#     if (opp_e.aos - opp_s.los) > MAX_SLEW_TIME:
#         return True

#     # Compute start and end line of sight vectors
#     z_start = _compute_los_vector(opp_s, opp_s.los)
#     z_end   = _compute_los_vector(opp_e, opp_e.aos)

#     # Compute slew time between orientations
#     t_slew = slew_time_single_axis(z_start, z_end)

#     # Check if slew time is less than the time available to complete the slew
#     if t_slew <= (opp_e.aos - opp_s.los):
#         return True
#     else:
#         return False

end # Constraints module