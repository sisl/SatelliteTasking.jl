__precompile__(true)
module MDP

# Julia Packages
using Random

using SatelliteDynamics.Time: Epoch
using SatelliteTasking.DataStructures: Image, Opportunity


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

# State constructor from image length
function MDPState(epc::Opportunity, image_list::Array{Image,1})    
    return MDPState(epc, BitArray(false for _ in 1:length(image_list)))
end

"""

"""
function mdp_compute_actions(state::MDPState, opportunities::Array{Opportunity, 1}, 
    constraint_list::Array{Function, 1}; breadth=0::Real)

    # NOTE: Assumes a time-ascending ordered order of opportunities

    # Get current opportunitity
    current_opp = state.time

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
function mdp_reward(state::MDPState, position_lookup::Dict{<:Integer, Image})
    
    r = 0.0
    
    for (i, img) in enumerate(state.locations)
        if img
            r += position_lookup[i].reward
        else
            r -= position_lookup[i].reward
        end
    end
    
    return r
end

function forward_search(state::MDPState, d::Integer, opportunities::Array{Opportunity, 1},
            constraint_list::Array{Function, 1},
            image_lookup::Dict{<:Integer, Image}, position_lookup::Dict{Image, <:Integer}; probabilities::Union{Dict{Opportunity, <:Real}, Nothing}=nothing,
            depth::Real=3, breadth::Real=3, gamma::Real=0.7)

    if d == 0
        return nothing, 0
    end

    # Initialize current optimal action and value
    astar, vstar = nothing, -Inf 

    for a in mdp_compute_actions(state, opportunities, constraint_list, breadth=breadth)
        v = mdp_reward(state, image_lookup)
        
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
            # Make action to take next image
            astar = future_states[1]

            # Simple vstar
            vstar = 0

            # More complex vstar
            success = MDPState(astar, deepcopy(state.locations))
            success.locations[position_lookup[astar.location]] = true
            vstar = mdp_reward(success, image_lookup)
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
    depth::Real=3, breadth::Real=3, gamma::Real=0.7)

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

end # End MDP module