# Export
export mdp_fs
export satellite_plan_mdp_fs

function mdp_depth_first_search(problem::PlanningProblem, state::MDPState, depth::Integer)
    
    # Early exit 
    if depth == 0
        return Noop(t_start=state.time), 0.0
    end

    astar, vstar = Noop(t_start=state.time), -Inf

    for a in mdp_state_actions(problem, state)
        v = mdp_reward(problem, state, a)

        for sp in mdp_reachable_states(problem, state, a)
            ap, vp =  mdp_depth_first_search(problem, sp, depth-1)

            # Strict Markov Update
            # v = v + problem.solve_gamma*vp

            # Get distance action is in the future
            t_diff = abs(ap.t_start - state.time)

            # Semi-markov update
            v = v + problem.solve_gamma^(t_diff)
        end

        if v > vstar
            astar, vstar = a, v
        end
    end

    return astar, vstar
end

function mdp_fs(problem::PlanningProblem, state; sat_id::Integer=1)

    # println("Foward Search Step")
    action, value = mdp_depth_first_search(problem, state, problem.solve_depth)

    if typeof(action) == Collect
        if action.location in state.requests
            println("Planned duplicate collect: $(action) - $(action.location)")
        end
    end

    # Step to next state with selected action
    state = mdp_step(problem, state, action)

    return state, action, value
end

function satellite_plan_mdp_fs(problem::PlanningProblem; sat_id::Integer=1)

    # Sort opportunities
    sort!(problem.opportunities, by = x -> x.t_start)
    
    # Set initial state
    init_opp = problem.opportunities[1]
    state = MDPState(time=init_opp.t_start, last_action=init_opp)

    plan = Opportunity[state.last_action]
    reward = 0.0

    while true
        state, action, value = mdp_fs(problem, state)

        # println("State: $state")
        # println("Action: $action")
        # println("Actions in horizon: $(length(find_actions(problem, state, problem.solve_horizon)))")
        # println("Value: $value")

        # Add state and action to plan
        push!(plan, action)
        reward += value

        # Break from search if terminal state
        if typeof(action) == Done
            break
        end
    end    

    return plan, reward
end