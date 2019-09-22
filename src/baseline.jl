# Exports
export satellite_plan_baseline
export satellite_plan_resource_baseline

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

function baseline_step(problem::SatPlanningProblem, state::SatMDPState, action::Opportunity)
    # If done return same state, just done
    if typeof(action) == Done
        return SatMDPState(time=state.time,
                last_action=Done(t_start=state.time),
                requests=copy(state.requests),
                power=state.power,
                data=state.data)
    end

    # Current state time
    time0 = state.time
    time  = state.time

    # Store last viable action
    last_action = state.last_action

    # Resource
    power = state.power
    data  = state.data

    # Copy observed requests
    requests = copy(state.requests)

    # Finish state
    done = false

    # Advance time to next action
    time = action.t_start

    # Update data generate
    power_generated = 0.0
    data_generated  = problem.spacecraft[1].datagen_backorbit*(time - time0)

    if typeof(action) == Noop
        # Do nothing if Noop
    elseif typeof(action) == Sunpoint
        # Charge for duration of sunpoint
        power_generated += problem.spacecraft[1].powergen_sunpoint*(time - time0)

    elseif typeof(action) == Collect
        # Update last action
        last_action = action

        collect_generation = action.duration * problem.spacecraft[1].datagen_image

        if data + collect_generation < 1.0
            # Only perform collect if we have capacilty
            data_generated  += collect_generation
            power_generated += action.duration * problem.spacecraft[1].powergen_image

            push!(requests, action.location)
        end

    elseif typeof(action) == Contact
        # Update last action
        last_action = action

        data_generated  += action.duration * problem.spacecraft[1].datagen_downlink
        power_generated += action.duration * problem.spacecraft[1].powergen_downlink
    else
        throw(ErrorException("Action $(string(action)) is not of a known type."))
    end

    # Update power and data numbers
    power += power_generated
    data  += data_generated

    # Ensure resource limits are enforced
    if power > 1.0 power = 1.0 end

    if power < 0.0 
        power = 0.0 
        done = true
    end

    if data > 1.0 data = 1.0 end
    if data < 0.0 data = 0.0 end

    return SatMDPState(time=time,
            last_action=last_action,
            requests=requests,
            power=power,
            data=data)
    
end

function next_feasible_action(problem::SatPlanningProblem, state::SatMDPState)
    
    # Find index of next possible action
    idx_start = find_next_opp_index(problem::SatPlanningProblem, state.time)

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

function satellite_plan_baseline(problem::SatPlanningProblem)

    # Get Initial State
    init_opp = problem.opportunities[1]
    state = SatMDPState(time=init_opp.t_start, last_action=init_opp)

    states = SatMDPState[state]
    plan = Opportunity[problem.opportunities[1]]
    requests = Request[]

    while next_feasible_action(problem, state) != nothing
        action = next_feasible_action(problem, state)

        # Advance state
        state = baseline_step(problem, state, action)

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

function satellite_plan_resource_baseline(problem::SatPlanningProblem; allow_repeats::Bool=false)

    # Get Initial State
    init_opp = problem.opportunities[1]
    state = SatMDPState(time=init_opp.t_start, last_action=init_opp)
    action = state.last_action

    states = SatMDPState[state]
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
        state = baseline_step(problem, state, action)

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