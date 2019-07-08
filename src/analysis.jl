# Exports
export analyze_plan

function analyze_plan(problem::PlanningProblem, plan::Array{<:Opportunity, 1})
    
    # Evaluate any constraint violations
    feasible = true
    for i in 1:length(plan)-1
        col_start = plan[i]
        col_end   = plan[i+1]

        # Only check validity of opportunity transitions
        if typeof(col_start) == Opportunity && typeof(col_end) == Opportunity
            # Set transition valid by default
            fe = true

            for constraint in problem.constraints
                # Use logical and to evaulate path feasibility on graph
                valid = valid && constraint(col_start, col_end)
            end

            if !valid
                println("Collect Failed:\n Start: $col_start\n End: $col_end\n")
                feasible = false
            end
        end
    end

    unique_requests = Request[]

    n_contacts = 0
    n_requests = 0
    n_dup_requests = 0

    reward = 0.0

    for opp in plan
        if typeof(opp) == Collect
            n_requests += 1
            reward += opp.reward
            if !(opp.location in unique_requests)
                push!(unique_requests, opp.location)
            else
                n_dup_requests += 1
            end
        elseif typeof(opp) == Contact
            n_contacts += 1
        end
    end
    n_unique_reqests = length(unique_requests)
    
    return feasible, reward, n_unique_reqests, n_dup_requests, n_requests, n_contacts
end