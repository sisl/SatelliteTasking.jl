# Exports
export get_abs_time, get_rel_time

"""
Computes the absolute time
"""
function get_abs_time(mdp::MDPProblem, t_since::Real)
    return mdp.t_start + t_since
end

"""
Computes the relative time
"""
function get_rel_time(mdp::MDPProblem, epc::Epoch)
    return epc - mdp.t_start
end