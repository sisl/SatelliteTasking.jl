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
function sp_construct_graph(collects::Array{Collect, 1}, constraint_list::Array{Function, 1}, horizon=0::Real)
    
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
            elseif horizon > 0 && end_collect.sow > (start_collect.sow + horizon)
                # Since we know collects are sorted we can break building the
                # transition grarph if the distance from the next to the next start
                # is greater than the look-ahead horizon
                break
            else
                # Set transition valid by default
                valid_transition = true

                for constraint in constraint_list
                    # Use logical and to evaulate path feasibility on graph
                    valid_transition = valid_transition && constraint(start_collect, end_collect)

                    # If not valid transition break early
                    if !valid_transition
                        break
                    end
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
function sp_solve_graph(graph::Dict{Collect, Array{Collect, 1}}, image_lookup::Dict{UUID, Image}, allow_repeats=false::Bool)
    optimal_path = Dict{Collect, Union{Union{Collect, Nothing}, Float64, Array{UUID}}}

    # Initialize nodes in optimal path
    for node in keys(graph)
        optimal_path[node] = (nothing, image_lookup[node.image_id].reward, [node.image_id])
    end
end


# def graph_compute_longest_path(graph, repeat_images=False):
#     optimal_path = {}

#     # Initialize Optimal Path
#     for node in graph.keys():
#         optimal_path[node] = (None, node.image.reward, [node.image])

#     # Iterate over nodes updating optimal weights
#     for node_i in optimal_path.keys():
#         for node_j in graph[node_i]:
#             if optimal_path[node_i][1] + node_j.image.reward > optimal_path[node_j][1]:
#                 if repeat_images:
#                     # If repeat images are allowed always update if reward is higher
#                     optimal_path[node_j] = (node_i, optimal_path[node_i][1] + node_j.image.reward, optimal_path[node_i][2] + [node_j.image])
#                 elif not repeat_images and node_j.image not in optimal_path[node_i][2]:
#                     # If repeat images are not allowed only update if the image is unique
#                     optimal_path[node_j] = (node_i, optimal_path[node_i][1] + node_j.image.reward, optimal_path[node_i][2] + [node_j.image])

#     # Backtrack to calculate optimal path (sort the nodes in the optimal path by highest reward)
#     end_node = sorted(zip(optimal_path.keys(), optimal_path.values()), key=lambda x: x[1][1], reverse=True)[0]

#     reward     = end_node[1][1] # Reward from optimal path
#     image_list = end_node[1][2] # Images observed on optimal path
#     path       = [end_node[0]]  # Optimal Path (terminating node)

#     while path[-1] != None:
#         path.append(optimal_path[path[-1]][0])

#     # Delete the earliest point in the list (since the second-to-last already points to the path start)
#     del path[-1]

#     # Reverse the path and image list
#     path.reverse()
#     image_list.reverse()

#     return path, reward, image_list

export sp_graph_policy
"""

"""
function sp_graph_policy(collects::Array{Collect, 1}, constraint_list; horizon=0::Real, allow_repeats=false::Bool)
    # Construct graph
    graph  = sp_construct_graph(collects, constraint_list, horizon)

    # Solve graph for taskign policy
    policy = sp_solve_graph(graph, allow_repeats=allow_repeats)

    return policy
end


end # End Graph module