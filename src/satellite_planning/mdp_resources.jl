__precompile__(true)
module MDPResource

# Julia Packages
using Random
using Printf

using SatelliteDynamics.Time: Epoch
using SatelliteTasking.DataStructures: Image, Opportunity


export MDPResourceState
"""
MDPResourceState 

Attributes:
- `time::Epoch` Current planning time
- `locations::BitArray{1}` Boolean bit array of whether a specific image has 
    been vistited. Index indicates image.
"""
mutable struct MDPResourceState
    time::Opportunity
    power::Float64
    data::Float64
    locations::BitArray{1}
end

# State constructor from image length
function MDPResourceState(epc::Opportunity, power::Real, data::Real, image_list::Array{Image,1})    
    return MDPResourceState(epc, BitArray(false for _ in 1:length(image_list)), power, data)
end

function Base.show(io::IO, state::MDPResourceState)

    seen = ""

    for l in state.locations
        if l == true
            seen = seen * "1"
        else
            seen = seen * "0"
        end
    end

    s = @sprintf "MDPResourceState(Opportunity: %s, Power: %.2f, Data: %.2f Locations: %s)" string(UInt64(pointer_from_objref(state.time)), base=16) state.power state.data seen

    print(io, s)
end

# Comparison operators for MDP state, relies on same allocation of opportunities
function Base.:(==)(state_left::MDPResourceState, state_right::MDPResourceState)
    return ((state_left.time == state_right.time) && 
            (state_left.power == state_right.power) &&
            (state_left.locations == state_right.locations))
end

function Base.:(!=)(state_left::MDPResourceState, state_right::MDPResourceState)
    return !(state_left == state_right)
end

##################
# Forward Search #
##################

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
function reachable_states(state::MDPResourceState, action::Opportunity, position_lookup::Dict{Image, <:Integer})
    
    
    # Failure to take collect would be same response at an updated time
    failure      = MDPResourceState(action, state.power, state.data, deepcopy(state.locations))
    
    # Success would be transitioning the reward function.
    success = MDPResourceState(action, failure.power, failure.data, deepcopy(failure.locations))
    success.locations[position_lookup[action.location]] = true
    
    return MDPResourceState[success, failure]
end

"""
Compute the reward for being in a given state.
"""
function mdp_reward(state::MDPResourceState, image_lookup::Dict{<:Integer, Image}, position_lookup::Dict{Image, <:Integer})
    
    r = 0.0
    
    for (i, img) in enumerate(state.locations)
        if img
            r += image_lookup[i].reward
        else
            r -= image_lookup[i].reward
        end
    end

    if state.power < 0
        r += -100000 # Large negative reward
    end

    if state.data > 1.0
        r += -10000 # Large negative reward for filling up storage
    end
    
    return r
    # return state.time.location.reward*state.locations[position_lookup[state.time.location]]
end

function forward_search(state::MDPResourceState, d::Integer, opportunities::Array{Opportunity, 1},
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
                    success = MDPResourceState(astar, deepcopy(state.locations))
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

function mdp_forward_step(state::MDPResourceState, action::Opportunity,
            position_lookup::Dict{Image, <:Integer};
            probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing)
    
    # Initialize next state
    state_next = MDPResourceState(action, state.power, state.data, deepcopy(state.locations))

    # Compute power generated over step time
    ts = state.time.eow
    te = state_next.time.sow

    # orb = state.orbit.
    
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
    state    = MDPResourceState(init_opp, 1.0, 0.0, images)

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

end # End MDP module