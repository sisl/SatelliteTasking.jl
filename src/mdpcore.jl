
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

export precompute_action_space
function precompute_action_space(problem::SatPlanningProblem)

    # First compute feasible transition lookup
    lt_feasible_transitions = Dict{Tuple{Integer, Integer}, Bool}()

    println("Computing Feasibility Lookup Table:")
    for so_idx in 1:length(problem.opportunities)
        so_opp = problem.lt_opportunities[so_idx]

        if so_idx % 100 == 0
            println("Current Index: $so_idx")
        end

        for eo_idx in (so_idx+1):length(problem.opportunities)
            eo_opp = problem.lt_opportunities[eo_idx]

            # Check feasibility constraints
            valid = true

            for constraint in problem.constraints
                if valid == false
                    break
                end

                valid = valid && constraint(problem.lt_opportunities[so_idx], problem.lt_opportunities[eo_idx])
            end

            # Push 
            lt_feasible_transitions[(so_idx, eo_idx)] = valid
        end
    end

    println("Current Index: $(length(problem.opportunities))")

    # Dict mapping (last collect/downlink action) -> {all possible actions}
    sunpoint_actions = Tuple{Opportunity, Integer}[]
    lt_feasible_actions = Dict{Tuple{Integer, Integer}, AbstractVector{Opportunity}}()

    # Below might work but it doesn't solve the problem of having sunpointed actions
    # be part of the 
    println("Computing Action Space Lookup:")
    for ca_idx in 1:length(problem.opportunities)
        ca_opp = problem.lt_opportunities[ca_idx]

        if ca_idx % 100 == 0
            println("Current Index: $ca_idx")
        end

        # Create sunpoint action once for time
        if ca_idx < length(problem.opportunities)
            sp_action = Sunpoint(id=length(problem.opportunities)+length(sunpoint_actions)+1, t_start=problem.opportunities[ca_idx+1].t_start)
            push!(sunpoint_actions, (sp_action, ca_idx+1)) # Record equivalent index of sunpoint action
        end

        for la_idx in 1:ca_idx-1
            la_opp = problem.lt_opportunities[la_idx]

            # Create Array for (current action, last action)
            lt_feasible_actions[(la_opp, ca_opp)] = Opportunity[]

            for fa_idx in (ca_idx+1):length(problem.opportunities)
                if lt_feasible_transitions[(la_idx, fa_idx)] == true
                    push!(lt_feasible_actions[(la_opp.id, ca_opp.id)], problem.opportunities[fa_idx])
                end

                if problem.solve_breadth > 0 && length(lt_feasible_actions[(la_opp.id, ca_opp.id)]) > problem.solve_breadth
                    break
                end
            end

            # Add sunpoint action as feasible action
            if ca_idx < length(problem.opportunities)
                push!(lt_feasible_actions[(la_opp.id, ca_opp.id)], sp_action)
            end
        end
    end

    println("Current Index: $(length(problem.opportunities))")

    println("Computing Sunpoint Indices:")
    for (sp_action, ca_idx) in sunpoint_actions
        ca_opp = problem.lt_opportunities[ca_idx]

        if ca_idx % 100 == 0
            println("Current Sunpoint Index: $ca_idx")
        end
        
        for la_idx in 1:ca_idx-1
            la_opp = problem.lt_opportunities[la_idx]

            # Copy Exisiting feasible actions 
            lt_feasible_actions[(la_opp.id, sp_action.id)] = copy(lt_feasible_actions[(la_opp.id, ca_opp.id)]) 
        end
    end

    println("Current Index: $(length(problem.opportunities))")

    # Update problem 
    problem.lt_feasible_actions = lt_feasible_actions
    problem.actions = vcat(problem.opportunities, [x[1] for x in sunpoint_actions])
    sort!(problem.actions, by = x -> x.t_start)
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

# State reward function
function POMDPs.reward(problem::SatPlanningProblem, state::SatMDPState, action::Opportunity)
    reward = 0.0

    # Update data generate
    power_generated = 0.0
    data_generated  = problem.spacecraft[1].datagen_backorbit*(action.t_start - state.time)

    if typeof(action) == Collect
        # Resources consumed by action
        data_generated  += action.duration * problem.spacecraft[1].datagen_image
        power_generated += action.duration * problem.spacecraft[1].powergen_image
        
        # Only reward if we have enough spare data space
        if (state.data + data_generated < 1.0)
            if !(action.location in state.requests)
                reward += action.reward
            end
            # elseif problem.solve_allow_repeats == true
            #     reward += action.reward*0.25
            # end
        else
            # Penalize collection when full up
            # reward -= action.reward
        end

    elseif typeof(action) == Contact
        # reward += 0.01*action.duration

        # Update data generation
        data_generated  += action.duration * problem.spacecraft[1].datagen_downlink
        power_generated += action.duration * problem.spacecraft[1].powergen_downlink
    elseif typeof(action) == Noop
        reward += 0.0
    elseif typeof(action) == Sunpoint
        # Charge for duration of sunpoint
        power_generated += problem.spacecraft[1].powergen_sunpoint*(action.t_start - state.time)

        # Have tiny bit of reward for charging
        # reward += (action.t_start - state.time)/state.last_action.duration
        reward += 0.0001*(action.t_start - state.time)
    end

    # Update resources to see if violations occur
    power = state.power + power_generated
    data  = state.data  + data_generated

    # Penalize running out of power
    if power <= 0
        reward -= 10000
    end

    # Don't penalize overcharge
    # if power > 1
    #     reward -= 10.0
    # end

    # Don't penalize full capacity
    # if data > 1
    #     reward -= 10.0
    # end

    # Don't penalize having everything down
    # if data < 0
    #     reward -= 10.0
    # end

    return reward
end

# Action Space
# POMDPs.actions(problem::SatPlanningProblem) = problem.actions
# function POMDPs.actions(problem::SatPlanningProblem, state::SatMDPState)
#     return mdp.lt_feasible_actions[(state.last_action.id, state.current_action.id)]
# end

# # Initial state funciton
# function POMDPs.initialstate(problem::SatPlanningProblem, rng::AbstractRNG)

#     # Get Initial Opportunity
#     init_opp = problem.opportunities[1]
    
#     # Create Initial State
#     state = SatMDPState(time=init_opp.t_start, last_action=init_opp)

#     return state
# end