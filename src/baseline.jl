# Exports
export satellite_plan_baseline
export satellite_plan_resource_baseline

function next_feasible_action(problem::PlanningProblem, state::MDPState)
    
    # Find index of next possible action
    idx_start = bs_find_firstgte(problem::PlanningProblem, state.time)

    if idx_start > length(problem.opportunities)
        return nothing
    end

    # Ensure index is for next action after current state
    while problem.opportunities[idx_start].t_start <= state.time
        idx_start += 1

        if idx_start > length(problem.opportunities)
            return nothing
        end
    end
    
    # Return all opportunities inside planning horizon
    fopps = problem.opportunities[idx_start:end]

    # List of all possible actions
    actions = Opportunity[]

    # Populate list of possible actions
    for opp in fopps
        # Skip Collects until contact if high data
        if typeof(opp) == Collect && state.data >= 0.5
            continue
        end

        # Valuate transition to see if valid
        valid = true

        # Check to see if each transition constraint is valid
        for constraint in problem.constraints
            if valid == false
                break
            end

            valid = valid && constraint(state.last_action, opp)
        end

        # If valid transition add to edges
        if valid == true 
            return opp
        end
    end
    
    return nothing
end

function satellite_plan_baseline(problem::PlanningProblem)

    # Get Initial State
    init_opp = problem.opportunities[1]
    state = MDPState(time=init_opp.t_start, last_action=init_opp)

    states = MDPState[state]
    plan = Opportunity[problem.opportunities[1]]
    requests = Request[]

    while next_feasible_action(problem, state) != nothing
        action = next_feasible_action(problem, state)

        # Advance state
        state = mdp_step(problem, state, action)

        push!(states, state)
        push!(plan, action)
    end

    reward = 0.0

    for action in plan
        if typeof(action) == Collect
            reward += action.reward
        end
    end

    return plan, reward
end

function satellite_plan_resource_baseline(problem::PlanningProblem; allow_repeats::Bool=false)

    # Get Initial State
    init_opp = problem.opportunities[1]
    state = MDPState(time=init_opp.t_start, last_action=init_opp)
    action = state.last_action

    states = MDPState[state]
    plan = Opportunity[problem.opportunities[1]]
    requests = Request[]

    while next_feasible_action(problem, state) != nothing
        # Get next feasible action
        action = next_feasible_action(problem, state)
        
        # If able to take collect or downlink do it
        power_generated = 0.0
        data_generated  = problem.spacecraft[1].datagen_backorbit*(action.t_start - state.time)

        if typeof(action) == Collect
            # Resources consumed by action
            data_generated  += action.duration * problem.spacecraft[1].datagen_image
            power_generated += action.duration * problem.spacecraft[1].powergen_image
        elseif typeof(action) == Contact
            # Update data generation
            data_generated  += action.duration * problem.spacecraft[1].datagen_downlink
            power_generated += action.duration * problem.spacecraft[1].powergen_downlink
        end
        
        if 0.0 < state.power + power_generated && state.data + data_generated < 1.0
            # Take action if able (e.g. do nothing)
        else
            # If not able sunpoint
            action = Sunpoint(t_start=action.t_start)
        end

        # println("State: $state")
        # println("Action: $action")

        # Advance state
        state = mdp_step(problem, state, action)

        push!(states, state)
        push!(plan, action)
    end

    reward = 0.0

    for action in plan
        if typeof(action) == Collect
            reward += action.reward
        end
    end

    return states, plan, reward
end