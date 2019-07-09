# Exports
export MDPState
export mdp_fs
export mdp_reward
export satellite_plan_mdp_fs

@with_kw struct MDPState
    time::Union{Epoch, Float64}
    last_action::Opportunity
    requests::Array{Request, 1} = Request[]
    power::Float64 = 1.0
    data::Float64 = 0.0
    done::Bool = false
end

function Base.hash(state::MDPState, h::UInt)
    return hash(state.time, hash(length(state.images), hash(state.power, hash(state.data, hash(state.done, hash(:Epoch, h))))))
end

function mdp_reward(problem::PlanningProblem, state::MDPState, action::Opportunity)
    reward = 0.0

    if typeof(action) == Collect && !(action.location in state.requests)
        # Resources consumed by action
        data_generated  = action.duration*problem.spacecraft[1].datagen_image
        power_generated = action.duration*problem.spacecraft[1].powergen_image

        # Only reward if we haven't overfilled
        if (state.data + data_generated < 1.0)
            reward += action.reward*(1.0 + problem.reward_alpha)^(-abs(action.t_start - state.time))
        end
    elseif typeof(action) == Contact
        reward += 0.01*action.duration
    end

    if state.power <= 0
        reward -= 1000000
    end

    return reward
end

function mdp_step(problem::PlanningProblem, state::MDPState, action::Opportunity)
    # If done return same state, just done
    if action == Done
        return MDPState(time=state.time,
                last_action=state.last_action,
                requests=copy(state.requests),
                power=state.power,
                data=state.data,
                done=true)
    end

    # Current state time
    time0 = state.time
    time  = state.time

    # Resource
    power = state.power
    data  = state.data

    # Copy observed requests
    requests = copy(state.requests)

    # Finish state
    done = false

    if typeof(action) == Noop
        # Advance time to next action without changing resources
        time = action.t_start
    elseif typeof(action) == Sunpoint
        time = action.t_start
    elseif typeof(action) == Collect
        time = action.t_start

        # Only perform collect if we have capacilty
        data_generated  = action.duration*problem.spacecraft[1].datagen_image
        power_generated = action.duration*problem.spacecraft[1].powergen_image

        if data + data_generated < 1.0
            power += power_generated
            data  += data_generated

            push!(requests, action.location)
        end

    elseif typeof(action) == Contact
        data_generated  = action.duration*problem.spacecraft[1].datagen_downlink
        power_generated = action.duration*problem.spacecraft[1].powergen_downlink
    else
        throw(ErrorException("Action $(string(action)) is not of a known type."))
    end

    # Ensure resource limits are enforced
    if power > 1.0 power = 1.0 end
    if power < 0.0 power = 0.0 end

    if data > 1.0 data = 1.0 end
    if data < 0.0 data = 0.0 end

    # If no other possible future actions problem state is done 
    if length(filter(x -> x.t_start > state.time, problem.opportunities)) == 0
        done = true
    end

    return MDPState(time=time,
            last_action=state.last_action,
            requests=requests,
            power=state.power,
            data=state.data,
            done=done)
    
end

function mdp_state_actions(problem::PlanningProblem, state::MDPState)
    # Limit search to future opportunities
    fopps = filter(x -> x.t_start > state.time, problem.opportunities)

    if length(fopps) == 0
        return Opportunity[Done]
    end

    # List of all possible actions
    actions = Opportunity[]

    # Populate list of possible actions
    for opp in fopps
        if problem.action_breadth > 0 && length(actions) == problem.action_breadth
            break
        end

        if typeof(state.last_action) == Collect
            # Valuate transition to see if valid
            valid = true

            # Check to see if each transition constraint is valid
            for constraint in problem.constraints
                if valid == false
                    break
                end

                valid = valid && constraint(start_opp, end_opp)
            end

            # If valid transition add to edges
            if valid == true
                push!(actions, opp)
            end
        else
            push!(actions, opp)
        end
    end

    # Add default actions at the end not to mess width action breadth-limit
    push!(actions, Noop(t_start=fopps[1].t_start))
    push!(actions, Sunpoint(t_start=fopps[1].t_start))

    # Return set of possible actions
    return actions
end

function mdp_reachable_states(problem::PlanningProblem, state::MDPState, action::Opportunity)
    return MDPState[mdp_step(problem, state, action)]
end

function mdp_depth_first_search(problem::PlanningProblem, state::MDPState, depth::Integer)
    
    # Early exit 
    if depth == 0
        return Noop, 0.0
    end

    astar, vstar = Noop, -Inf

    for a in mdp_state_actions(problem, state)
        v = mdp_reward(problem, state, a)

        for sp in mdp_reachable_states(problem, state, a)
            ap, vp = mdp_depth_first_search(problem, state, depth-1)

            v = v + problem.solve_gamma*vp
        end

        if v > vstar
            astar, vstar = a, v
        end
    end

    return astar, vstar
end

function mdp_fs(problem::PlanningProblem, state; sat_id::Integer=1)

    action, value = mdp_depth_first_search(problem, state, problem.solve_depth)

    # Step to next state with selected action
    state = mdp_step(problem, state, action)
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

        println("State: $state")
        println("Action: $action")
        println("Value: $value")

        # Add state and action to plan
        push!(plan, state.last_action)
        reward += value

        # Break from search if terminal state
        if state.done == true
            break
        end
    end    

    return plan, reward
end