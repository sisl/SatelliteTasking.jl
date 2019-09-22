
###########################
# Action Space Definition #
###########################

"""
Binary search method to find the first opportunity with an index greater than
the target simulation time (Epoch or elapsed time)
"""
function find_next_opp_index(problem::SatPlanningProblem, target::Union{Epoch, Real})
    # Set search indices
    l = 1
    r = length(problem.opportunities)
    m = 0
    
    # Perform binary search to get first opportunity with start time
    # greater or equal to the target value
    while l < r
        m = floor(Int, (l+r)/2)
        # println("l: $l, r: $r, avg: $((l+r)/2), m: $m, t: $(problem.opportunities[m].t_start)")
        
        if problem.opportunities[m].t_start < target
            l = m + 1
        else
            r = m
        end
    end
    
    return m
end

"""
Find all possible (not necessarily feasible) actions given a state and a
planning horizon.
"""
function find_actions(problem::SatPlanningProblem, epc::Epoch, horizon::Real)
    # Find index of next possible action
    idx_start = find_next_opp_index(problem::SatPlanningProblem, epc)

    # Ensure index is for next action after current state
    while problem.opportunities[idx_start].t_start <= epc
        idx_start += 1
    end

    # Find index of first action outside of horizon
    idx_end = find_next_opp_index(problem::SatPlanningProblem, epc+horizon)
    
    # Get index of last action inside horizon
    idx_end = idx_end - 1
    
    # Return all opportunities inside planning horizon
    return problem.opportunities[idx_start:idx_end]
end

function upcoming_actions(problem::SatPlanningProblem, epc::Epoch)
    # Limit search to future opportunities
    fopps = find_actions(problem, epc, problem.solve_horizon)

    # List of all possible actions
    actions = Opportunity[]

    # Populate list of possible actions
    for opp in fopps
        if problem.solve_breadth > 0 && length(actions) == problem.solve_breadth
            break
        end

        push!(actions, opp)
    end

    # Add default actions at the end not to mess width action breadth-limit
    # Set time to earliest time of next action
    if length(fopps) > 0
        push!(actions, Sunpoint(t_start=fopps[1].t_start))
        # push!(actions, Noop(t_start=next_action_time))
    end

    # Return set of possible actions
    return actions
end

export precompute_actions
function precompute_actions(problem::SatPlanningProblem)

    # First compute feasible transition lookup
    lt_feasible_transitions = Dict{Tuple{Opportunity, Opportunity}, Bool}()

    for so_idx in 1:length(problem.opportunities)
        for eo_idx in (so_idx+1):length(problem.opportunities)

        end
    end


    # Dict mapping (last collect/downlink action) -> {all possible actions}
    possible_actions = Dict{Opportunity, AbstractVector{Opportunity}}()

    sunpoint_actions = Tuple{Opportunity, Integer}[]
    feasible_actions = Dict{Tuple{Opportunity, Opportunity}, AbstractVector{Opportunity}}()


    # Below might work but it doesn't solve the problem of having sunpointed actions
    # be part of the 
    for ca_idx in 1:length(problem.opportunities)
        ca_opp = problem.lt_opportunities[ca_idx]

        if ca_idx % 100 == 0
            println("Current Index: $ca_idx")
        end

        for la_idx in 1:ca_idx-1
            la_opp = problem.lt_opportunities[la_idx]
            # println("s$ca_idx, $la_idx")

            # Create Array for (current action, last action)
            feasible_actions[(ca_opp, la_opp)] = Opportunity[]

            for fa_idx in (ca_idx+1):length(problem.opportunities)
                # Check feasibility constraints
                valid = true

                for constraint in problem.constraints
                    if valid == false
                        break
                    end

                    valid = valid && constraint(problem.lt_opportunities[la_idx], problem.lt_opportunities[fa_idx])
                end

                if valid == true
                    push!(feasible_actions[(ca_opp, la_opp)], problem.lt_opportunities[fa_idx])
                end
            end

            # Create Sunpoint Opportunities and copy feasible actions
            if ca_idx < length(problem.opportunities)
                # Create action
                sp_action = Sunpoint(t_start=problem.opportunities[ca_idx+1].t_start)

                push!(sunpoint_actions, (sp_action, ca_idx+1)) # Record equivalent index of sunpoint action
                push!(feasible_actions[(ca_opp, la_opp)], sp_action)
            end
        end
    end

    println("Current Index: $(length(problem.opportunities))")

    for (sp_action, ca_idx) in sunpoint_actions
        if ca_idx % 100 == 0
            println("Current Sunpoint Index: $ca_idx")
        end

        for la_idx in 1:ca_idx-1
            la_opp = problem.lt_opportunities[la_idx]
            
            eq_opp = problem.lt_opportunities[ca_idx]

            # Copy Exisiting feasible actions 
            # push!(feasible_actions[(sp_action, la_opp)], feasible_actions[(eq_opp, la_opp)])
        end
    end

    println("Current Index: $(length(problem.opportunities))")

    # Update problem 
    # problem.lt_feasible_actions = feasible_actions
    problem.actions = vcat(problem.opportunities, sunpoint_actions)
    sort!(problem.actions, by = x -> x.t_start)

    return problem.actions

    # for opp in problem.opportunities
    #     # Initialize 
    #     possible_actions[opp] = Opportunity[]

    #     for idx in (opp.id+1):length(problem.opportunities)
    #         # Check feasibility constraints
    #         valid = true

    #         for constraint in problem.constraints
    #             if valid == false
    #                 break
    #             end

    #             valid = valid && constraint(opp, problem.lt_opportunities[idx])
    #         end

    #         if valid == true
    #             push!(possible_actions[opp], problem.lt_opportunities[idx])
    #         end
    #     end

    #     # Add Sunpoint and Noop actions
    #     if opp.id < length(problem.opportunities)
    #         # Create action
    #         sp_action = Sunpoint(t_start=problem.opportunities[opp.id+1].t_start)

    #         push!(possible_actions[opp], sp_action)
    #         # push!(possible_actions[opp], Noop(t_start=problem.opportunities[opp.id+1].t_start))
    #     end
    # end

    # n_opps = 0
    # n_actions = 0
    # for idx in 1:length(problem.opportunities)-1
    #     opp = problem.opportunities[idx]
    #     actions = upcoming_actions(problem, opp.t_start)
    #     println("$(opp.t_start) - $actions")
    #     n_opps += 1
    #     n_actions += 1
    # end

    # return possible_actions
end

export precompute_feasible_actions
function precompute_feasible_actions(problem::SatPlanningProblem)

end

###########################
# Action Space Definition #
###########################

# Terminal States
# POMDPs.isterminal(mdp::TaskingMDP, s::TaskingState) = (length(mdp.action_lookup[a.sow]) == 0)
POMDPs.isterminal(problem::SatPlanningProblem, s::SatMDPState) = s.done
POMDPs.discount(problem::SatPlanningProblem) = 1.0

# # Generate next state from state and action
# function POMDPs.generate_s(problem::SatPlanningProblem, s::SatMDPState, a::Opportunity, rng::AbstractRNG)
#     # Calculate Charging
#     t    = s.last_collect.eow
#     te   = a.sow
#     step = 1.0

#     # Power generated
#     pg = 0.0

#     while t < te
#         pg += eclipse_conical(t, state(tle, t))
#         t  += step
#     end
    
#     # New state values
#     time  = a.sow
#     last_collect = a
#     power = min(max(s.power + pg - mdp.power_rate_image*a.location.collect_duration, 0.0), 1.0)
#     data  = s.data # Do nothing with data for now
#     done  = s.done
#     dead  = s.dead
    
#     if power <= 0.0
#         dead = true
#         done = true # We're done if we're dead
#     end
    
#     # If no futuer actions exit
#     if length(mdp.action_lookup[a.sow]) == 0
#         done = true
#     end
    
#     # Copy image array and set taken to true
#     image_array = copy(s.image_array)
#     image_array[mdp.position_lookup[a.location]] = true
    
#     # Initialize next state
#     sp = SatMDPState(time, last_collect, power, data, image_array, dead, done)
    
#     return sp
# end

# # State reward function
# function POMDPs.reward(p::SatPlanningProblem, s::SatMDPState, a::Opportunity)
#     # Ensure satellite still has pwoer
#     if s.dead == true
#         # Big penalty
#         return -1000000.0
#     end
    
#     # Reward as sum of all current images
#     r = 0.0
#     for (i, img) in enumerate(s.image_array)
#         if img == true
#             r += p.image_lookup[i].reward
#         else
#             r -= p.image_lookup[i].reward
#         end
#     end
    
#     return r
# end

# # Action Space
# POMDPs.actions(problem::SatPlanningProblem) = mdp.opportunities
# function POMDPs.actions(problem::SatPlanningProblem, s::SatMDPState)
#     return mdp.action_lookup[s.time]
# end

# # Initial state funciton
# function POMDPs.initialstate(problem::SatPlanningProblem, rng::AbstractRNG)

#     # Get Initial Opportunity
#     init_opp = problem.opportunities[1]
    
#     # Create Initial State
#     state = SatMDPState(time=init_opp.t_start, last_action=init_opp)

#     return state
# end