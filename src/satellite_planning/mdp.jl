__precompile__(true)
module MDP

# Julia Packages
using Random
using Printf

using SatelliteDynamics.Time: Epoch
using SatelliteTasking.DataStructures: Image, Opportunity


export MDPState
"""
MDPState 

Attributes:
- `time::Epoch` Current planning time
- `locations::BitArray{1}` Boolean bit array of whether a specific image has 
    been vistited. Index indicates image.
"""
mutable struct MDPState
    time::Opportunity
    locations::BitArray{1}
end

function Base.show(io::IO, state::MDPState)

    seen = ""

    for l in state.locations
        if l == true
            seen = seen * "1"
        else
            seen = seen * "0"
        end
    end

    s = @sprintf "MDPState(Opportunity: %s, Locations: %s)" string(UInt64(pointer_from_objref(state.time)), base=16) seen

    print(io, s)
end

# Comparison operators for MDP state, relies on same allocation of opportunities
function Base.:(==)(state_left::MDPState, state_right::MDPState)
    return ((state_left.time == state_right.time) && (state_left.locations == state_right.locations))
end

function Base.:(!=)(state_left::MDPState, state_right::MDPState)
    return ((state_left.time != state_right.time) || (state_left.locations != state_right.locations))
end

# State constructor from image length
function MDPState(epc::Opportunity, image_list::Array{Image,1})    
    return MDPState(epc, BitArray(false for _ in 1:length(image_list)))
end

function mdp_compute_actions(current_opp::Opportunity, opportunities::Array{Opportunity, 1}, 
    constraint_list::Array{Function, 1}; breadth=0::Real)

    # NOTE: Assumes a time-ascending ordered order of opportunities

    # Limit search to future opportunities
    future_opportunities = filter(x -> x.sow > current_opp.eow, opportunities)

    # Populate list of possible actions
    actions = Opportunity[]

    for end_opportunity in future_opportunities
        # Exit early if we've found enough actions for our breadth search
        if breadth > 0 && length(actions) == breadth
            break
        end

        # Check if valid transition
        valid_transition = true
        for constraint in constraint_list
            # If not valid transititon break early
            if valid_transition == false
                continue
            end

            # Use logical AND ot evaluate continued feasibility
            valid_transition = valid_transition && constraint(current_opp, end_opportunity)
        end

        # If valid add action to list of valid state actions
        if valid_transition
            push!(actions, end_opportunity)
        end
    end

    return actions
end

export extract_images
"""
Given an array of opportunities extract an array containing all unique images
which are being observed.

Arguments:
    - `opportunities::Array{Opportunity, 1}`: All observation opportunities

Returns:
    - `images::Array{Image, 1}`: Array of unique images
"""
function extract_images(opportunities::Array{Opportunity, 1})
    images = Image[]

    # Iterate through opportunities adding unique images
    for opp in opportunities
        if !(opp.location in images)
            push!(images, opp.location)
        end
    end

    return images
end

"""
Given an array of locations createt a lookup dictionary that maps the image to an
index in a array.

Arguments:
    - `images::Array{Image, 1}`: Array of Image locations

Returns:
    - `image_lookup::Dict{Integer, Image}` Dictionary mapping array indices to 
        image objects.
    - `position_lookup::Dict{Integer, Image}` Dictionary mapping array indices to 
        image objects.
"""
function create_lookups(images::Array{Image, 1})
    image_lookup    = Dict(idx => image for (idx, image) in enumerate(images))
    position_lookup = Dict(image => idx for (idx, image) in enumerate(images))

    return image_lookup, position_lookup
end

"""
Enumerate possible states reachable 
"""
function reachable_states(state::MDPState, action::Opportunity, position_lookup::Dict{Image, <:Integer})
    # Failure to take collect would be same response at an updated time
    failure      = MDPState(state.time, deepcopy(state.locations))
    failure.time = action
    
    # Success would be transitioning the reward function.
    success = MDPState(failure.time, deepcopy(failure.locations))
    success.locations[position_lookup[action.location]] = true
    
    return MDPState[success, failure]
end

"""
Compute the reward for being in a given state.
"""
function mdp_reward(state::MDPState, image_lookup::Dict{<:Integer, Image}, position_lookup::Dict{Image, <:Integer})
    
    r = 0.0
    
    for (i, img) in enumerate(state.locations)
        if img
            r += image_lookup[i].reward
        else
            r -= image_lookup[i].reward
        end
    end
    
    return r
    # return state.time.location.reward*state.locations[position_lookup[state.time.location]]
end

function forward_search(state::MDPState, d::Integer, opportunities::Array{Opportunity, 1},
            constraint_list::Array{Function, 1},
            image_lookup::Dict{<:Integer, Image}, position_lookup::Dict{Image, <:Integer}; probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing,
            depth::Real=3, breadth::Real=3, gamma::Real=0.95)

    if d == 0
        return nothing, 0
    end

    # Initialize current optimal action and value
    astar, vstar = nothing, -Inf 

    for a in mdp_compute_actions(state.time, opportunities, constraint_list, breadth=breadth)
        v = mdp_reward(state, image_lookup, position_lookup)
        
        for sp in reachable_states(state, a, position_lookup)
            # Continue forward search of space
            ap, vp = forward_search(sp, d-1, opportunities, constraint_list, 
                        image_lookup, position_lookup, probabilities=probabilities,
                        depth=depth, breadth=breadth, gamma=gamma)

            # Update state reward. Probabilities is the transition function directly
            if probabilities != nothing
                v = v + gamma*probabilities[a]*vp
            else
                v = v + gamma*vp
            end
        end

        if v > vstar
            astar, vstar = a, v
        end
    end

    # It is possible there are no actions found. If a reachable state still
    # exists take the next one
    if astar == nothing
        future_states = filter(x -> x.sow > state.time.eow, opportunities)
        if length(future_states) > 0
            # Look through possible future transitions and take the first feasible 
            for fstate in future_states
                # Check if valid transition
                valid_transition = true
                for constraint in constraint_list
                    # If not valid transititon break early
                    if valid_transition == false
                        continue
                    end

                    # Use logical AND ot evaluate continued feasibility
                    valid_transition = valid_transition && constraint(state.time, fstate)
                end
                
                if valid_transition
                    # Make action to take next image
                    astar = future_states[1]

                    # Simple vstar
                    vstar = 0

                    # More complex vstar
                    success = MDPState(astar, deepcopy(state.locations))
                    success.locations[position_lookup[astar.location]] = true
                    vstar = mdp_reward(success, image_lookup, position_lookup)

                    # Successful transition found, break
                    break
                end
            end
        end
    end

    return astar, vstar
end

function mdp_forward_step(state::MDPState, action::Opportunity,
            position_lookup::Dict{Image, <:Integer};
            probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing)
    
    # Initialize next state
    state_next = MDPState(action, deepcopy(state.locations))

    
    if probabilities != nothing
        # Tranition state probabilistically if probabilities are available
        if probabilities[action] <= rand()
            state_next.locations[position_lookup[action.location]] = true
        end
    else    
        # Set state to opportunityed
        state_next.locations[position_lookup[action.location]] = true
    end

    return state_next
end

# MDP forward search algorithm
export mdp_solve_forward_search
"""
Solve MDP using basic forward search algorithm.
"""
function mdp_solve_forward_search(opportunities::Array{Opportunity, 1}, 
    constraint_list::Array{Function, 1}, probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing;
    depth::Real=3, breadth::Real=3, gamma::Real=0.95)

    # Extract image list
    images = extract_images(opportunities)

    # Compute image lookup
    image_lookup, position_lookup = create_lookups(images)

    # Set Initial state
    init_opp = opportunities[collect(keys(opportunities))[findmin(collect([o.sow for o in opportunities]))[2]]]
    state    = MDPState(init_opp, images)

    # Store plan
    plan = Union{Opportunity, Nothing}[]

    # Take initial step
    action, value = forward_search(state, depth, opportunities, 
                        constraint_list, image_lookup, position_lookup, 
                        probabilities=probabilities, breadth=breadth, gamma=gamma)

    # Add action to plan
    push!(plan, action)

    # Transition state forward
    state = mdp_forward_step(state, action, position_lookup)
    # state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)

    while true
        action, value = forward_search(state, depth, opportunities, 
                        constraint_list, image_lookup, position_lookup, 
                        probabilities=probabilities, breadth=breadth, gamma=gamma)

        if action == nothing
            # Exit early
            break
        end

        # Add action to plan
        push!(plan, action)

        # Transition state forward
        state = mdp_forward_step(state, action, position_lookup)
        # state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)
    end

    # Compute plan reward
    reward = 0
    image_list = Image[]
    
    for (i, image_collected) in enumerate(state.locations)
        if image_collected
            image   = image_lookup[i]
            reward += image.reward
            push!(image_list, image)
        end
    end
        
    return plan, reward, image_list
end

###########################
# Monte Carlo Tree Search #
###########################

function mcts_generate_rollout(opportunities::Array{Opportunity, 1}, constraint_list::Array{Function, 1}; breadth::Integer=10)

    rollout_policy = Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}}()

    for start_opportunity in opportunities

        # Initialize Array
        rollout_policy[start_opportunity] = Array{Tuple{Opportunity, <:Real}, 1}[]

        # Populate list of possible actions
        actions = mdp_compute_actions(start_opportunity, opportunities, constraint_list, breadth=breadth)

        num_actions = length(actions)
        for a in actions
            push!(rollout_policy[start_opportunity], (a, 1.0/num_actions))
        end
    end

    return rollout_policy
end

function mcts_get_optimal_action(state::MDPState, N::Dict{Opportunity, Dict{Opportunity, Int32}}, Q::Dict{Opportunity, Dict{Opportunity, Float64}}; c=1.0)
    # Find optimal action
    amax, vmax = nothing, -Inf

    Ns = sum([N[state.time][as] for as in keys(Q[state.time])])
    for a in keys(Q[state.time])

        if Ns == 0
            # If Ns is zero use this otherwise log(Ns) = -Inf
            v = Q[state.time][a]
        else
            v = Q[state.time][a] + c*sqrt(log(Ns)/N[state.time][a])
        end

        if v > vmax
            amax, vmax = a, v
        end
    end

    return amax
end

function mcts_tree_search(state::MDPState, d::Integer, rollout_policy::Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}};
            opportunities::Array{Opportunity, 1},
            probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing,
            constraint_list::Array{Function, 1}, 
            position_lookup::Dict{Image, <:Integer},
            image_lookup::Dict{<:Integer, Image}, 
            T::Array{Opportunity, 1},
            N::Dict{Opportunity, Dict{Opportunity, Int32}},
            Q::Dict{Opportunity, Dict{Opportunity, Float64}},
            c::Real=1.0, gamma::Real=0.95, breadth::Integer=10)
   
    if d == 0
        return 0.0
    end

    # If state not seen, add it and rollout policy 
    if !(state.time in T)
        # println("Encountered new state $(string(UInt64(pointer_from_objref(state.time)), base=16)).")

        # Initiallies N & Q arrays
        N[state.time] = Dict{Opportunity, Int32}()
        Q[state.time] = Dict{Opportunity, Float64}()

        for a in mdp_compute_actions(state.time, opportunities, constraint_list, breadth=breadth)
            # Initialize state action reward value if not
            N[state.time][a] = 0
            Q[state.time][a] = 0
        end

        # Add state to observed states
        push!(T, state.time)

        # Return reward for state using rollout_policy
        return mcts_rollout(state, d, rollout_policy, 
                gamma=gamma, 
                probabilities=probabilities, 
                position_lookup=position_lookup, 
                image_lookup=image_lookup)
    else    
        # println("Encountered state $(string(UInt64(pointer_from_objref(state.time)), base=16)) - $(sum([N[state.time][as] for as in keys(Q[state.time])])) times previously.")
    end

    # Find optimal action
    amax = mcts_get_optimal_action(state, N, Q, c=c)    

    # No future actions return nothing
    if amax == nothing
        return 0.0
    end

    # Simulate next state
    next_state = mdp_forward_step(state, amax, position_lookup, probabilities=probabilities)
       
    # Reward for action
    # r = amax.location.reward*next_state.locations[position_lookup[amax.location]]
    r = mdp_reward(next_state, image_lookup, position_lookup)

    q = r + gamma*mcts_tree_search(next_state, d-1, rollout_policy, 
                    opportunities=opportunities, probabilities=probabilities,
                    constraint_list=constraint_list,
                    position_lookup=position_lookup,
                    image_lookup=image_lookup, 
                    T=T, Q=Q, N=N, c=c, gamma=gamma, breadth=breadth)
    
    N[state.time][amax] = N[state.time][amax] + 1
    Q[state.time][amax] = Q[state.time][amax] + (q - Q[state.time][amax])/N[state.time][amax]

    return q
end

function sample_rollout(state::MDPState, rollout_policy::Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}})
    # Generate random number
    p = rand()

    # Set cumulative probability
    cum_prob = 0.0

    for ap in rollout_policy[state.time]
        cum_prob += ap[2]

        if p <= cum_prob
            return ap[1]
        end
    end

    return nothing
end

function mcts_rollout(state, d, rollout_policy::Dict{Opportunity, Array{Tuple{Opportunity, <:Real}}}; 
            probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing,
            position_lookup::Dict{Image, <:Integer},
            image_lookup::Dict{<:Integer, Image}, 
            gamma::Real=0.95)
            
    if d == 0
        return 0
    end

    # Sample rollout policy
    action = sample_rollout(state, rollout_policy)

    # If no possible future states return 0
    if action == nothing
        return 0
    end

    # Simulate next state
    next_state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)
    
    # Reward for action
    # r = action.location.reward*next_state.locations[position_lookup[action.location]]
    r = mdp_reward(next_state, image_lookup, position_lookup)

    return r + gamma*mcts_rollout(next_state, d-1, rollout_policy, 
                        gamma=gamma,
                        probabilities=probabilities, 
                        position_lookup=position_lookup,
                        image_lookup=image_lookup)
end


export mdp_solve_mcts
"""
Solve MDP witht Monte Carlo Tree Search
"""
function mdp_solve_mcts(opportunities::Array{Opportunity, 1}, 
    constraint_list::Array{Function, 1}, probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing;
    depth::Real=10, breadth::Integer=10, gamma::Real=0.95, c::Real=0.75, max_iterations::Integer=10)

    # Extract image list
    images = extract_images(opportunities)

    # Compute image lookup
    image_lookup, position_lookup = create_lookups(images)

    # Initialize MCTS variables
    plan = Union{Opportunity, Nothing}[]
    T = Opportunity[]                                    # Visited States 
    N = Dict{Opportunity, Dict{Opportunity, Int32}}()    # State-Action Exploration Count
    Q = Dict{Opportunity, Dict{Opportunity, Float64}}()  # State-Action Reward

    # Generate rollout policy
    rollout_policy = mcts_generate_rollout(opportunities, constraint_list, breadth=breadth)

    # Set Initial state
    init_opp = opportunities[collect(keys(opportunities))[findmin(collect([o.sow for o in opportunities]))[2]]]
    state    = MDPState(init_opp, images)

    # println("Current state: $(state.time.sow)")

    # Run MCTS for fix number of iterations
    for i in 1:max_iterations
        # println("Tree search Iteration $i")
        q = mcts_tree_search(state, depth, rollout_policy, 
                    opportunities=opportunities,
                    probabilities=probabilities,
                    constraint_list=constraint_list,
                    position_lookup=position_lookup,
                    image_lookup=image_lookup,
                    T=T,
                    Q=Q,
                    N=N,
                    c=c,
                    gamma=gamma,
                    breadth=breadth)
    end

    # Get Optimal action for state
    action = mcts_get_optimal_action(state, N, Q, c=0)

    # Add action to plan
    push!(plan, action)

    # Transition state forward probabilistically
    state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)
    # println("Current state: $(state.time.sow)")

    while true
        # Run MCTS for fix number of iterations
        for i in 1:max_iterations
            # println("Tree search Iteration $i")
            q = mcts_tree_search(state, depth, rollout_policy, 
                        opportunities=opportunities,
                        probabilities=probabilities,
                        constraint_list=constraint_list,
                        position_lookup=position_lookup,
                        image_lookup=image_lookup,
                        T=T,
                        Q=Q,
                        N=N,
                        c=c,
                        gamma=gamma,
                        breadth=breadth)
        end

        # Get Optimal action for state
        action = mcts_get_optimal_action(state, N, Q, c=0)

        if action == nothing
            @debug "No next action at time: $(state.time) Stopping search. $(mdp_compute_actions(state.time, opportunities, constraint_list, breadth=breadth))"
            # Exit early
            break
        end

        # Add action to plan
        push!(plan, action)

        # Transition state forward probabilistically
        state = mdp_forward_step(state, action, position_lookup, probabilities=probabilities)
        @debug "Current state: $(state.time.sow)"
    end


    # Compute plan reward
    reward = 0
    image_list = Image[]
    
    for (i, image_collected) in enumerate(state.locations)
        if image_collected
            image   = image_lookup[i]
            reward += image.reward
            push!(image_list, image)
        end
    end
        
    return plan, reward, image_list
end

end # End MDP module