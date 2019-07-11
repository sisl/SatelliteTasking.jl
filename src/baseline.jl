# Exports
export satellite_plan_baseline

function next_feasible_action(problem::PlanningProblem, last_action::Opportunity, requests::Array{Request, 1}; allow_repeats::Bool=false)
    
    future_opps = filter(x -> x.t_start > last_action.t_start, problem.opportunities)

    for opp in future_opps
        # Skip repeats disallowed
        if allow_repeats == false && typeof(opp) == Collect && opp.location in requests
            continue
        end

        valid = true

        # Check to see if each transition constraint is valid
        for constraint in problem.constraints
            if valid == false
                break
            end

            valid = valid && constraint(last_action, opp)
        end

        # If valid transition add to edges
        if valid == true 
            return opp
        end
    end
    
    return nothing
end

function satellite_plan_baseline(problem::PlanningProblem; allow_repeats::Bool=false)

    plan = Opportunity[problem.opportunities[1]]
    requests = Request[]

    while next_feasible_action(problem, plan[end], requests, allow_repeats=allow_repeats) != nothing
        push!(plan, next_feasible_action(problem, plan[end], requests, allow_repeats=allow_repeats))
        push!(requests, plan[end].location)
    end

    reward = 0.0

    for action in plan
        if typeof(action) == Collect
            reward += action.reward
        end
    end

    return plan, reward
end