# Exports
export get_abs_time, get_rel_time

"""
Computes the absolute time
"""
function get_abs_time(problem::SatPlanningProblem, t_since::Real)
    return problem.t_start + t_since
end

"""
Computes the relative time
"""
function get_rel_time(problem::SatPlanningProblem, epc::Epoch)
    return epc - problem.t_start
end