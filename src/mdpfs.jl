# Export
export mdp_fs
export satellite_plan_mdp_fs

function mdp_reachable_states(problem::SatPlanningProblem, state::SatMDPState, action::Opportunity)
    return SatMDPState[baseline_step(problem, state, action)]
end

function mdp_depth_first_search(problem::SatPlanningProblem, state::SatMDPState, depth::Integer)
    
    # Early exit 
    if depth == 0
        return Noop(t_start=state.time), 0.0
    end

    astar, vstar = Done(t_start=state.time), -Inf

    for a in problem.lt_feasible_actions[(state.last_cdo_action.id, state.last_action.id)]
        v = POMDPs.reward(problem, state, a)

        for sp in mdp_reachable_states(problem, state, a)
            ap, vp = mdp_depth_first_search(problem, sp, depth-1)

            # Strict Markov Update
            # v = v + problem.solve_discount*vp

            # Get distance action is in the future
            t_diff = abs(ap.t_start - state.time)

            # Semi-markov update
            v = v + vp*problem.solve_discount^(t_diff)
        end

        if v > vstar
            astar, vstar = a, v
        end
    end

    return astar, vstar
end

function mdp_fs(problem::SatPlanningProblem, state::SatMDPState)

    # println("Foward Search Step")
    action, value = mdp_depth_first_search(problem, state, problem.solve_depth)

    if typeof(action) == Collect
        if action.location.id in state.request_ids
            # println("Planned duplicate collect: $(action) - $(action.location)")
        end
    end

    # if typeof(action) == Contact
    #     println("Taking Contact: $action")
    # end

    if typeof(action) != Done && action.t_start == state.time
        # throw(ErrorException("Took action that isn't advancing time...\n $state - $action"))
        @debug "Took action that isn't advancing time...\n $state - $action"
        return state, Done(t_start=state.time), 0
    end

    # Step to next state with selected action
    state = baseline_step(problem, state, action)

    return state, action, value
end

function satellite_plan_mdp_fs(problem::SatPlanningProblem)

    # Sort opportunities
    sort!(problem.opportunities, by = x -> x.t_start)
    
    # Set initial state
    init_opp = problem.opportunities[1]
    state = SatMDPState(time=init_opp.t_start, last_action=init_opp, last_cdo_action=init_opp)

    states = SatMDPState[state]
    plan = Opportunity[state.last_action]
    total_reward = 0.0

    while true
        state, action, reward = mdp_fs(problem, state)

        # Break from search if terminal state
        if typeof(action) == Done
            break
        end

        # Add state and action to plan
        push!(states, state)
        push!(plan, action)
        total_reward += reward
    end    

    return states, plan, total_reward
end