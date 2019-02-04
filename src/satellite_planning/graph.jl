__precompile__(true)
module Graph

using UUIDs

# Package Imports
using SatelliteTasking.DataStructures: Image, Opportunity, Collect

export sp_construct_graph
"""

Arguments:
- `collects::Array{Collect,1}` Array of collects to plan the optimal tasking schedule for
- `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
- `horizon::Real` Look-ahead horizon for constructing graphs. Transition further than this apart are not considered. Not used if 0

Returns:
- `collect_graph` Graph of feasible collect transitions
"""
function sp_construct_graph(collects::Array{Collect, 1}, constraint_list::Array{Function, 1}; horizon=0::Real)
    
    # Sort Collects to ensure they are in time-asecnding order
    sort!(collects, by = x -> x.sow)

    # Initialize validate transition graph
    graph = Dict{Collect, Array{Collect, 1}}()

    for start_collect in collects
        # List of valid edges/transitions for start_collect
        graph[start_collect] = Array{Collect, 1}[]
        
        for end_collect in collects
            # Decide to insert edge
            if (start_collect == end_collect || 
                start_collect.image == end_collect.image ||
                end_collect.sow < start_collect.eow)
                # Skip insertion if same collect, image, or starts before the current
                # collection ends (no instantaneous manevers)
                continue
            elseif horizon > 0 && end_collect.sow > (start_collect.eow + horizon)
                # Since we know collects are sorted we can break building the
                # transition grarph if the distance from the next to the next start
                # is greater than the look-ahead horizon
                break
            else
                # Set transition valid by default
                valid_transition = true

                for constraint in constraint_list
                    # If not valid transition break early
                    if valid_transition == false
                        continue
                    end

                    # Use logical and to evaulate path feasibility on graph
                    valid_transition = valid_transition && constraint(start_collect, end_collect)
                end

                if valid_transition
                    push!(graph[start_collect], end_collect)
                end
            end
        end
    end
    
    return graph
end

# Solve graph 
export sp_solve_graph
"""
Solve graph for decision policy

Arguments:
- `graph::Dict{Collect, Array{Collect, 1}}`
- `allow_repeats::Bool` Allow images to be collected multiple times over the course of a plan

Returns:
- `collect_policy::Array{Collect}` List of collects to take in the order which they should be taken
"""
function sp_solve_graph(graph::Dict{Collect, Array{Collect, 1}}; allow_repeats=false::Bool)
    optimal_path = Dict{Collect, Tuple{Union{Collect, Nothing}, Float64, Array{Image}}}()

    # Initialize nodes in optimal path
    for node in keys(graph)
        optimal_path[node] = (nothing, node.image.reward, [node.image])
    end

    # Iterate over nodes updating optimal weights
    for node_i in sort!(collect(keys(optimal_path)), by = x -> x.sow) # Subtle, but this needs to be sorted in time order to give optimal policy
        for node_j in sort!(graph[node_i], by = x -> x.sow) # Also needs to be sorted for same reason, and Julia dicts don't order keys by default
            if (optimal_path[node_i][2] + node_j.image.reward) > optimal_path[node_j][2]
                if allow_repeats == true
                    # If repeat images are allowed always update if reward is higher
                    optimal_path[node_j] = (node_i, optimal_path[node_i][2] + node_j.image.reward, push!(copy(optimal_path[node_i][3]), node_j.image))
                elseif allow_repeats == false && !(node_j.image in optimal_path[node_i][3])
                    # If repeat images are not allowed only update if the image is unique
                    optimal_path[node_j] = (node_i, optimal_path[node_i][2] + node_j.image.reward, push!(copy(optimal_path[node_i][3]), node_j.image))
                end
            end
        end
    end
    
    # Find node with highest reward:
    current_max = 0.0
    end_node = nothing
    for node in keys(graph)
        if optimal_path[node][2] > current_max
            current_max = optimal_path[node][2]
            end_node    = node
        end
    end

    @debug end_node

    # end_node = sort!(collect(zip(keys(optimal_path), values(optimal_path))), by = x -> x[2][2], rev=true)[1]

    # Extract values
    reward     = optimal_path[end_node][2] # Reward from optimal path
    image_list = optimal_path[end_node][3] # Images observed on optimal path
    path       = Union{Collect, Nothing}[end_node]  # Optimal Path (terminating node)

    while path[end] != nothing
        push!(path, optimal_path[path[end]][1])
    end
    
    # Remove last element since it's going to be nothing
    pop!(path)

    return (path, reward, image_list)
end

export sp_graph_policy
"""
Solve for optimal collect plan using a graph traversal algorithm.

Arguments:
- `collects::Array{Collect,1}` Array of collects to plan the optimal tasking schedule for
- `constraints::Array{Any, 1}` Array of constraint function which may limit feasible transitions
- `horizon::Real` Look-ahead horizon for constructing graphs. Transition further than this apart are not considered. Not used if 0
- `allow_repeats::Bool` Allow images to be collected multiple times over the course of a plan

Returns:
- `collect_policy::Array{Collect}` List of collects to take in the order which they should be taken
"""
function sp_graph_policy(collects::Array{Collect, 1}, constraint_list; horizon=0::Real, allow_repeats=false::Bool)
    # Construct graph
    graph  = sp_construct_graph(collects, constraint_list, horizon=horizon)

    # Solve graph for taskign policy
    path, reward, image_list = sp_solve_graph(graph, allow_repeats=allow_repeats)

    return path, reward, image_list
end


end # End Graph module