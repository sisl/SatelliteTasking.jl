__precompile__(true)
module MDP

# Julia Packages
using Random
using Printf

using SatelliteDynamics.Time: Epoch, mjd
using SatelliteDynamics.OrbitDynamics: eclipse_conical
using SatelliteDynamics.SGPModels: TLE
using SatelliteTasking.DataStructures: Image, GroundStation, Opportunity, Orbit


########################
# General MDP Planning #
########################

# Default spacecraft model
sc_model_default = Dict{String,Float64}(
    "power_generation" => 0.0,
    "power_draw_downlink" => 0.0,
    "power_draw_imaging"  => 0.0,
    "data_generation_backorbit" => 0.0,
    "data_generation_imaging" => 0.0,
    "data_downlink" => 0.0
)

export MDPState
"""
MDPState 

Attributes:
- `time::Epoch` Current planning time
- `locations::BitArray{1}` Boolean bit array of whether a specific image has 
    been vistited. Index indicates image.
"""
mutable struct MDPState
    time::Epoch
    last_action::Union{Opportunity, Symbol}
    images::Array{Image, 1}
    downlink_queue::Vector{Image}
    power::Float64
    data::Float64
    done::Bool
end

function Base.show(io::IO, state::MDPState)

    s = @sprintf "MDPState(time: %s, power: %.2f, data: %.2f, image_count: %d)" state.time state.power state.data length(state.images)

    print(io, s)
end

# Comparison operators for MDP state, relies on same allocation of opportunities
function Base.:(==)(state_left::MDPState, state_right::MDPState)
    return ((state_left.time == state_right.time) &&
            (state_left.last_action == state_right.last_action) &&
            (state_left.downlink_queue == state_right.downlink_queue) &&
            (state_left.images == state_right.images) &&
            (state_left.power == state_right.power) &&
            (state_left.data == state_right.data) &&
            (state_left.done == state_right.done))
end

function Base.:(!=)(state_left::MDPState, state_right::MDPState)
    return ((state_left.time != state_right.time) ||
            (state_left.last_action != state_right.last_action) ||
            (state_left.downlink_queue != state_right.downlink_queue) ||
            (state_left.images != state_right.images) ||
            (state_left.power != state_right.power) ||
            (state_left.data != state_right.data) ||
            (state_left.done != state_right.done))
end

# State constructor from image length
function mdp_compute_actions(state::MDPState, opportunities::Array{Opportunity, 1}; 
            constraint_list::Array{Function, 1}=Function[], breadth=0::Real)

    # Limit search to future opportunities
    future_opportunities = filter(x -> x.sow > state.time, opportunities)

    # Populate list of possible actions
    actions = Union{Opportunity, Symbol}[]

    # Add image collection opportunities
    for end_opportunity in future_opportunities
        # Exit early if we've found enough actions for our breadth search
        if breadth > 0 && length(actions) == breadth
            break
        end

        # Check if valid transition
        valid_transition = true
        if typeof(state.last_action) == Opportunity
            # Check constraints 
            for constraint in constraint_list
                # If not valid transititon break early
                if valid_transition == false
                    continue
                end

                # Use logical AND ot evaluate continued feasibility
                valid_transition = valid_transition && constraint(state.last_action, end_opportunity)
            end

            # If valid add action to list of valid state actions
            if valid_transition
                push!(actions, end_opportunity)
            end
        else
            # Last action was a Symbol
            if state.last_action == :SUNPOINT
                push!(actions, end_opportunity)
            end
        end
    end

    # Only add forced actions if other actions are left
    if length(future_opportunities) > 0
        # Add sunpoint action as valid action
        push!(actions, :SUNPOINT)
    end

    return actions
end

"""
Compute the reward for being in a given state.
"""
function mdp_reward(state::MDPState, action::Union{Symbol, Opportunity}; 
            alpha::Real=0.0, sc_model::Dict{String, Float64}=sc_model_default)

    # Bias function
    r = 0.0

    # println("R1: $r")
    if action == :SUNPOINT
        r += 1000
        # println("Rewarding SUNPOINT")
    end
    
    if typeof(action) == Opportunity && typeof(action.location) == Image
        # Only reward new images
        if !(action.location in state.images)
            # Resources consumed by action
            dg  = action.duration*sc_model["data_downlink"] 
            pg = action.location.collect_duration*sc_model["power_draw_imaging"]

            # Only reward if we haven't over-filled data
            if ((state.data  + dg) < 1.0) && ((state.power + pg) > 0.0)
                # println("Collecting: $(state.power) + $pg")
                r += action.location.reward*(1.0 + alpha)^(-abs(action.sow - state.time))
            else
                r += -10000
            end
        end
    elseif typeof(action) == Opportunity && typeof(action.location) == GroundStation
        # Reward for downlinking data
        r += 0.01*action.duration
    end

    # println("R2: $r")


    # Penalize running out tof power 
    if state.power <= 0.0
        r += -1000000
    end

    # println("R3: $r")

    # # Penalize 
    # if action == :SUNPOINT
    #     r += 1000
    #     println("Rewarding SUNPOINT")
    # end

    # println("State: $(state.time) - A: $action - Reward: $r")

    return r
end

function extract_images(opportunities::Array{Opportunity, 1})
    images = Image[]

    for opp in opportunities
        if !(opp.location in images) && typeof(opp.location) == Image
            push!(images, opp.location)
        end
    end

    return images
end

################################
# Deterministic Forward Search #
################################

function reachable_states(state::MDPState, action::Union{Symbol, Opportunity}, 
            opportunities::Array{Opportunity, 1},
            sc_model::Dict{String, Float64}=sc_model_default,
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}())
    # Compute all possible reachable states given acition
    # (In this case it's deterministic) 

    time0 = state.time
    time  = state.time

    # Spacecraft resources
    power = state.power
    data  = state.data

    # Shallow copy observed images
    images = copy(state.images)

    # Shallow copy downlink Queue
    dlqueue = copy(state.downlink_queue)

    if action == :SUNPOINT
        # Sunpoint transition
        # idx  = findfirst(o -> o.sow > time, opportunities)
        # time = opportunities[idx].sow
        time = next_action_time[time0]

        # Generate power for duration of step
        power += sc_model["power_generation"]*(time - time0)

        # Recharge
        # println("Investigating sunpoint state.")
        # println("Power before: $power")
        # power += sc_model["power_generation"]*(opportunities[idx].sow-state.time)
        # println("Power after: $power")
    elseif typeof(action) == Opportunity && typeof(action.location) == Image
        # Image Collection

        # Advance time
        time = action.sow

        # Only perform collect if we have data capacity
        data_generated = action.location.collect_duration*sc_model["data_generation_imaging"]

        # Check we don't overflow data buffer
        if (data + data_generated) < 1.0
            # Decrement power
            power += action.location.collect_duration*sc_model["power_draw_imaging"]

            # Increment data 
            data += data_generated

            push!(images, action.location)
        else
            # Decrement Power
            power += action.location.collect_duration*sc_model["power_draw_imaging"]

            # Set data to max
            data = 1.0
        end

    elseif typeof(action) == Opportunity && typeof(action.location) == GroundStation
        # Ground Contact

        # Advance time
        time = action.sow

        # Decrement downlink power rate
        power += 0

        # Decrement data        
        data += action.duration*sc_model["data_downlink"]

    end

    # Ensure Power limits
    if power > 1.0
        power = 1.0
    end
    if power < 0.0
        power = 0.0
    end

    # Ensure Data Limits
    if data > 1.0
        data = 1.0
    end
    if data < 0.0
        data = 0.0
    end

    sp = MDPState(time, action, images, dlqueue, power, data, false)

    return MDPState[sp]
end

function forward_search(state::MDPState, d::Int, 
            opportunities::Array{Opportunity, 1},
            sc_model::Dict{String, Float64}=sc_model_default; 
            constraint_list::Array{Function, 1}=Function[],
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}(),
            depth::Real=3, breadth::Real=3, gamma::Real=1.0, alpha::Real=0.0)

    if d == 0
        return :SUNPOINT, 0.0
    end

    # Initialize current optimal action and value
    astar, vstar = :DONE, -Inf 

    for a in mdp_compute_actions(state, opportunities, constraint_list=constraint_list, breadth=breadth)
        v = mdp_reward(state, a, alpha=alpha)
        # println("d: $d - a: $a - v: $v")
        
        for sp in reachable_states(state, a, opportunities, sc_model, next_action_time)
            # Continue forward search of space
            ap, vp = forward_search(sp, d-1, opportunities, sc_model,
                        constraint_list=constraint_list,
                        next_action_time=next_action_time,
                        depth=depth, breadth=breadth, gamma=gamma)
            
            # Update state reward.
            v = v + gamma*vp
        end

        if v > vstar
            astar, vstar = a, v
        end
    end

    return astar, vstar
end

function mdp_forward_step(state::MDPState, action::Union{Opportunity, Symbol}, 
            opportunities::Array{Opportunity, 1},
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}(),
            sc_model::Dict{String, Float64}=sc_model_default)
    
    # Current state time
    time0 = state.time
    time  = state.time

    # Resources
    power = state.power
    data  = state.data

    # Agent state
    done = false

    # Shallow copy observed images
    images = copy(state.images)

    # Shallow copy downlink Queue
    dlqueue = copy(state.downlink_queue)

    # Transition based on action type
    if action == :SUNPOINT
        # Sunpoint transition
        # idx  = findfirst(o -> o.sow > time, opportunities)
        # time = opportunities[idx].sow
        time = next_action_time[time0]
        duration = (time - time0)

        # Generate power for duration of step
        power += sc_model["power_generation"]*duration

    elseif typeof(action) == Opportunity && typeof(action.location) == Image
        # Opportunity Transition

        # Advance time
        time = action.sow

        # Only perform collect if we have data capacity
        data_generated = action.location.collect_duration*sc_model["data_generation_imaging"]

        # Check we don't overflow data buffer
        if (data + data_generated) < 1.0
            # Decrement power
            power += action.location.collect_duration*sc_model["power_draw_imaging"]

            # Increment data 
            data += data_generated

            push!(images, action.location)
        else
            # Decrement Power
            power += action.location.collect_duration*sc_model["power_draw_imaging"]

            # Set data to max
            data = 1.0
        end

    elseif typeof(action) == Opportunity && typeof(action.location) == GroundStation
        # Ground Contact

        # Advance time
        time = action.sow

        # Decrement downlink power rate
        power += 0

        # Decrement data        
        data += action.duration*sc_model["data_downlink"]

    else
        throw(ErrorException("Unknown action type $(string(action))"))
    end

    # Ensure Power limits
    if power > 1.0
        power = 1.0
    end
    if power < 0.0
        power = 0.0
    end

    # Ensure Data Limits
    if data > 1.0
        data = 1.0
    end
    if data < 0.0
        data = 0.0
    end

    # If no other possible future actions thte problem is done
    if length(filter(x -> x.sow > state.time, opportunities)) == 0
        done = true
    end

    # Initialize next state
    state_next = MDPState(time, action, images, dlqueue, power, data, done)

    return state_next
end

# MDP forward search algorithm
export mdp_solve_forward_search
"""
Solve MDP using basic forward search algorithm.
"""
function mdp_solve_forward_search(opportunities::Array{Opportunity, 1}, 
            constraint_list::Array{Function, 1}=Function[];
            depth::Real=3, breadth::Real=3, gamma::Real=1.0, alpha::Real=0.0, 
            sc_model::Dict{String, Float64}=sc_model_default)

    # Sort Opportunities
    sort!(opportunities, by = x -> x.sow)
    
    # Extract image list
    images = extract_images(opportunities)

    # Next Action list 
    next_action_time = Dict{Epoch, Epoch}()
    for idx in 1:(length(opportunities)-1)
        next_action_time[opportunities[idx].sow] = opportunities[idx+1].sow
    end

    # Set Initial state
    init_opp = opportunities[collect(keys(opportunities))[findmin(collect([o.sow for o in opportunities]))[2]]]
    state    = MDPState(init_opp.sow, init_opp, Image[init_opp.location], Vector{Image}(), 1.0, 0.0, false)

    # Store plan (action history) and state history
    plan   = Union{Opportunity, Symbol, Nothing}[]
    states = Union{MDPState, Nothing}[]

    # Take initial step
    action, value = forward_search(state, depth, opportunities, sc_model, 
                        constraint_list=constraint_list, depth=depth, 
                        breadth=breadth, gamma=gamma, alpha=alpha, next_action_time=next_action_time)

    # Add initial state to search
    push!(states, state)
    push!(plan, action)

    # Transition state forward
    state = mdp_forward_step(state, action, opportunities, next_action_time, sc_model)

    # println("State: $(string(state))")

    while true
        # Compute next action
        action, value = forward_search(state, depth, opportunities, sc_model, 
                            constraint_list=constraint_list, depth=depth, 
                            breadth=breadth, gamma=gamma, alpha=alpha, next_action_time=next_action_time)

        # println("State: $(state.time) - Action: $action - Reward: $value")
        # println("Next action: $(string(action))")

        if action == :DONE
            state.done = true
            break
        end

        # Add state and action to plan 
        push!(states, state)
        push!(plan, action)

        # Transition state forward
        state = mdp_forward_step(state, action, opportunities, next_action_time, sc_model)

        # println("State: $(string(state))")
    end

    # Compute plan reward
    image_list = copy(state.images)
    reward     = 0
    for image in image_list
        reward += image.reward
    end
        
    return states, plan, reward, image_list
end

####################
# Branch and Bound #
####################

###########################
# Monte Carlo Tree Search #
###########################

function policy_even_select(state::MDPState, opportunities::Array{Opportunity, 1};
            constraint_list::Array{Function, 1}=Function[], breadth=10::Real)
    
    # Compute all possible future actions
    actions = mdp_compute_actions(state, opportunities, constraint_list=constraint_list, breadth=breadth)

    if length(actions) == 0
        return :DONE
    end

    # Randomly sample action out of all possible actions
    # Generate random number
    p = rand()

    # Map propability back to action index uniformly
    action_idx = ceil(Int, rand()*length(actions)) 
    action = actions[action_idx]

    return action
end

function mcts_rollout(state::MDPState, d::Int, opportunities::Array{Opportunity, 1};
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}(),
            constraint_list::Array{Function, 1}=Function[],
            sc_model::Dict{String, Float64}=sc_model_default, 
            breadth::Int=10, gamma::Real=1.0, alpha::Real=1.0)
    # Return if no depth
    if d == 0
        return 0
    end

    # Sample rollout policty
    action = policy_even_select(state, opportunities, constraint_list=constraint_list, breadth=breadth)

    if action == :DONE
        return 0
    end

    # Calculate reward for being in state: state and taking action: action
    reward = mdp_reward(state, action, alpha=alpha, sc_model=sc_model)

    # Simulate next state
    next_state = mdp_forward_step(state, action, opportunities, next_action_time, sc_model)

    return reward + gamma*mcts_rollout(next_state, d-1, opportunities,
                        next_action_time=next_action_time,
                        constraint_list=constraint_list,
                        sc_model=sc_model,
                        breadth=breadth,
                        gamma=gamma,
                        alpha=alpha)
end

function mcts_simulate(state::MDPState, d::Integer,
            opportunities::Array{Opportunity, 1},
            T::Array{MDPState, 1},
            L::Dict{MDPState, Array{Union{Opportunity, Symbol}, 1}},
            N::Dict{Tuple{MDPState, Union{Opportunity, Symbol}}, Int32},
            Q::Dict{Tuple{MDPState, Union{Opportunity, Symbol}}, Float64};
            constraint_list::Array{Function, 1},
            sc_model::Dict{String, Float64}=sc_model_default,
            next_action_time::Dict{Epoch, Epoch}=Dict{Epoch, Epoch}(),
            breadth::Integer=10,
            c::Real=1.0, gamma::Real=1.0, alpha::Real=1.0)
   
    if d == 0
        return 0.0
    end

    # If state not seen, add it and rollout policy 
    if !(state in T)
        # println("Encountered new state $(string(UInt64(pointer_from_objref(state.time)), base=16)).")
        # println("Encounted new state: $state")
        
        # Add state to observed states
        push!(T, state)
        
        # Initialize state lookup if first observation
        L[state] = Union{Opportunity, Symbol}[]

        for action in mdp_compute_actions(state, opportunities, constraint_list=constraint_list, breadth=breadth)
            # if action == :SUNPOINT
            #     # println("Considering $(string(action)) action.")
            # end

            if !(action in L[state])
                push!(L[state], action)
            end

            # Initiallies N & Q arrays
            N[(state, action)] = 0
            Q[(state, action)] = 0.0
        end

        # Return reward for state using rollout_policy
        return mcts_rollout(state, d, opportunities,
                    next_action_time=next_action_time, 
                    sc_model=sc_model, breadth=breadth, gamma=gamma, alpha=alpha)
    end

    # Calculate N(s) by summing over all observed states 
    Ns = 0
    for a in L[state]
        Ns = Ns + N[(state, a)]
    end

    println("Ns: $Ns")

    # Get Action
    astar, vstar = :DONE, -Inf
    for a in L[state]
        println("S: $state - A: $a")
        println("Q: $(Q[(state, a)]) + c*sqrt(log(Ns)/N[(state, a)]: $(c*sqrt(log(Ns)/N[(state, a)]))")
        if (Q[(state, a)] + c*sqrt(log(Ns)/N[(state, a)])) > vstar
            astar, vstar = a, Q[(state, a)]
        end
    end

    # println("astar: $(typeof(astar))")

    # (randomly) sample next state and calculate reward for that action
    next_state = mdp_forward_step(state, astar, opportunities, next_action_time, sc_model)
    r = mdp_reward(state, astar, alpha=alpha, sc_model=sc_model)

    # Apply MCTS to get reward of next state
    q = r + gamma*mcts_simulate(state, d-1, opportunities, T, L, N, Q,
            constraint_list=constraint_list, sc_model=sc_model, next_action_time=next_action_time,
            breadth=breadth, c=c, gamma=gamma, alpha=alpha)

    # Update MCTS parameters 
    N[(state, astar)] = N[(state, astar)] + 1
    Q[(state, astar)] = Q[(state, astar)] + (q-Q[(state, astar)])/N[(state, astar)]

    return q
end

function mcts_optimal_action(state::MDPState,
    L::Dict{MDPState, Array{Union{Opportunity, Symbol}, 1}},
    Q::Dict{Tuple{MDPState, Union{Opportunity, Symbol}}, Float64})

    astar, vstar = :DONE, -Inf
    for a in L[state]
        # println("State: $(state.time) - Action: $a - Reward: $(Q[(state, a)])")
        if Q[(state, a)] > vstar
            astar, vstar = a, Q[(state, a)]
        end
    end

    # println("Optitmal action: $action")

    return astar
end

# MDP forward search algorithm
export mdp_solve_mcts
"""
Solve MDP using basic forward search algorithm.
"""
function mdp_solve_mcts(opportunities::Array{Opportunity, 1}, 
            constraint_list::Array{Function, 1}=Function[];
            depth::Real=3, breadth::Real=3, gamma::Real=1.0, c::Real=1.0, 
            iterations::Int=25,
            sc_model::Dict{String, Float64}=sc_model_default)

    # Sort Opportunities
    sort!(opportunities, by = x -> x.sow)
    
    # Extract image list
    images = extract_images(opportunities)

    # Next Action list 
    next_action_time = Dict{Epoch, Epoch}()
    for idx in 1:(length(opportunities)-1)
        # println("Added key $(opportunities[idx].sow)")
        next_action_time[opportunities[idx].sow] = opportunities[idx+1].sow
    end

    # Orbit
    orbit = opportunities[1].orbit

    # Set Initial state
    init_opp = opportunities[collect(keys(opportunities))[findmin(collect([o.sow for o in opportunities]))[2]]]
    state    = MDPState(init_opp.sow, init_opp, Image[init_opp.location], Vector{Image}(), 1.0, 0.0, false)

    # Store plan (action history) and state history
    plan   = Union{Opportunity, Symbol, Nothing}[]
    states = Union{MDPState, Nothing}[]

    # Initialize MCTS variables
    T = MDPState[]                                                   # Visited States
    L = Dict{MDPState, Array{Union{Opportunity, Symbol}, 1}}()       # State-Action Lookup
    N = Dict{Tuple{MDPState, Union{Opportunity, Symbol}}, Int32}()   # State-Action Exploration Count
    Q = Dict{Tuple{MDPState, Union{Opportunity, Symbol}}, Float64}() # State-Action Reward

    # Take initial step
    # for i in 1:iterations
    value = mcts_simulate(state, depth, opportunities, T, L, N, Q, 
                        constraint_list=constraint_list,
                        sc_model=sc_model,
                        breadth=breadth, 
                        gamma=gamma, 
                        next_action_time=next_action_time)
    # end

    # Add initial state to search
    action = mcts_optimal_action(state, L, Q)
    # println("At state $(state.time) - optimal action: $(typeof(action)) - Value: $value")
    state  = mdp_forward_step(state, action, opportunities, next_action_time, sc_model)
    push!(states, state)
    push!(plan, action)

    

    # Transition state forward
    while true
        # for i in 1:iterations
        value = mcts_simulate(state, depth, opportunities, T, L, N, Q, 
                        constraint_list=constraint_list,
                        sc_model=sc_model,
                        breadth=breadth,
                        gamma=gamma, 
                        next_action_time=next_action_time)
        # end

        # Add initial state to search
        action = mcts_optimal_action(state, L, Q)
        if action == :DONE
            state.done = true
            break
        end
        # println("At state $(state.time) - optimal action: $(typeof(action)) - Value: $value")

        state  = mdp_forward_step(state, action, opportunities, next_action_time, sc_model)

        if state.done
            println("Found DONE flag in state. Terminating early.")
            break
        end

        push!(states, state)
        push!(plan, action)
    end

    # Compute plan reward
    image_list = copy(state.images)
    reward     = 0
    for image in image_list
        reward += image.reward
    end
        
    return states, plan, reward, image_list
end

end