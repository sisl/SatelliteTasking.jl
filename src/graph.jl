# Exports
export construct_graph
export solve_graph
export satellite_plan_graph

"""
Construct decision graph for solving network
"""
function construct_graph(problem::PlanningProblem; sat_id::Integer=1)
    # Sort opportunities to ensure time ascending order
    sort!(problem.opportunities, by = x -> x.t_start)

    # Get max slew time
    sc = problem.spacecraft[sat_id]

    graph = Dict{Opportunity, Array{<:Opportunity, 1}}()

    for start_opp in problem.opportunities
        # Create list of valid transitions for start_collect
        graph[start_opp] = Opportunity[]

        for end_opp in problem.opportunities
            
            if (start_opp == end_opp) || (start_opp.t_start >= end_opp.t_end)
               # Conditions to skip insertion
            else
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
                    push!(graph[start_opp], end_opp)
                end

            end
        end
    end

    return graph 
end

"""
Solve a graph to find optimal transition path
"""
function solve_graph(problem::PlanningProblem, graph::Dict{Opportunity, Array{<:Opportunity, 1}};
            allow_repeats::Bool=false, contact_reward::Real=0.0)
   
    # Create graph based on incoming edges for each node
    igraph = Dict{Opportunity, Array{<:Opportunity, 1}}()

    for opp in problem.opportunities
        igraph[opp] = Opportunity[]
    end

    # For each node in the original graph add the source node as a possible
    # predecessor to each outgoing edge it has
    for node in sort(collect(keys(graph)), by = x -> x.t_start)
        for edge in graph[node]
            push!(igraph[edge], node)
        end
    end

    # Initialize structure to store optimal decision at each node
    optimal_path = Dict{Opportunity, Tuple{Union{Opportunity, Nothing}, Float64, Array{Request, 1}}}()

    for opp in problem.opportunities
        # Initialize node weight
        if typeof(opp) == Collect
            if allow_repeats == false
                optimal_path[opp] = (nothing, opp.reward, Request[opp.location])
            else
                optimal_path[opp] = (nothing, opp.reward, Request[])
            end
        else
            optimal_path[opp] = (nothing, 0.0, Request[])
        end
    end

    # For each node iterate over all incoming nodes
    for nodej in sort(collect(keys(igraph)), by = x -> x.t_start)
        # Iterate over all possible incoming nodes to node j
        for nodei in sort(igraph[nodej], by = x -> x.t_start)
            # Compute reward for current action
            rij = optimal_path[nodei][2]
            if typeof(nodej) == Collect
                # Calculate reward for taking j after the node i path
                if allow_repeats == true || !(nodej.location in optimal_path[nodei][3])
                    rij += nodej.reward
                end
            elseif typeof(nodej) == Contact
                rij += contact_reward
            end

            # If incoming node i is higher than anything seen previously set it
            # as new predecessory
            if rij > optimal_path[nodej][2]
                request_history = copy(optimal_path[nodei][3])
                if allow_repeats == false && typeof(nodej.location) == Request
                    push!(request_history, nodej.location)
                end
                optimal_path[nodej] = (nodei, rij, request_history)
            end
        end
    end

    return optimal_path
end

"""
Schedule satellite activity using graph planning. 

Returns:
- `plan::Array{<:Opportunity, 1}`
- `reward::Float64`
- `requests::{<:Request, 1}`
"""
function satellite_plan_graph(problem::PlanningProblem; sat_id::Integer=1,
            allow_repeats::Bool=false, contact_reward::Real=0)
    # Construct graph to solve 
    graph = construct_graph(problem, sat_id=sat_id)

    # Solve graph 
    optimal_path = solve_graph(problem, graph, 
                        allow_repeats=allow_repeats,
                        contact_reward=contact_reward)

    # Extract prameters of optimal path
    path     = Opportunity[]
    reward   = 0.0

    # Find node with highest reward:
    current_max = -Inf
    end_node = nothing
    for node in keys(optimal_path)
        if optimal_path[node][2] > current_max
            current_max = optimal_path[node][2]
            end_node    = node
        end
    end
    reward = current_max

    # Backtrack from end node to get optimal path
    push!(path, end_node)
    current_node = end_node
    while optimal_path[current_node][1] != nothing
        previous_node = optimal_path[current_node][1]
        push!(path, previous_node)
        current_node = previous_node
    end

    # Reverse path so it's in chronological order
    path = reverse!(path)

    return path, reward
end