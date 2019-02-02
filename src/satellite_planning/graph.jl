__precompile__(true)
module Graph

# Package Imports
using SatelliteTasking.DataStructures: Collect

export sp_construct_graph
"""

Arguments:
- `collects::Array{Collect,1}`
- `constraints`
- `horizon`

Returns:
- `collect_graph` Graph of feasible transition graph
"""
function sp_construct_graph(collects::Array{Collect, 1}, constraints, horizon=0::Real)
    
    # Initialize validate transition graph
    graph = Dict{Collect, Array{Collect, 1}}()

    # for start_collect in collects
    #     # List of valid edges/transitions for start_collect
    #     edges = Array{Collect, 1}[]
        
    #     for end_collect in collects
    #         # Decide to insert edge
    #         if (start_collect == end_collect || 
    #             start_collect.image_id == end_collect.image_id ||
    #             end_collect.sow < start_collect.sow)
    #             # Skip insertion if same collect, image, or is earlier in time
    #             continue
    #         elseif horizon > 0 && end_collect.sow > 

    #         end
    #     end
    # end
    
    return graph
end

# # Iterate over all nodes in graph constructing edges
# for opp_start in opportunities:
#     # List to store node edges (which nodes they are going for)
#     graph[opp_start] = []

#     for opp_end in opportunities:
#         if opp_start == opp_end or opp_start.image == opp_end.image:
#             # Skip if the same opportunity or same image
#             continue
#         elif opp_end.aos < opp_start.los:
#             # Skip if the end is before the start
#             continue
#         elif horizon and opp_end.aos > (opp_start.los + horizon):
#             # Condition to exit early is only considering transitions within a certain horizon
#             break
#         else:
#             # Transition is default valid
#             valid_transition = True

#             for constraint in constraint_list:
#                 # Use logical and to evaulate path feasibility on graph
#                 valid_transition = valid_transition and constraint(opp_start, opp_end)

#             if valid_transition:
#                 graph[opp_start].append(opp_end)
#             end

end # End Graph module