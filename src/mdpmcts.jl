# Exports
export MCTSState
export satellite_plan_mdp_mcts

@with_kw mutable struct MCTSState
    T = MDPState[]                                                   # Visited States
    L = Dict{MDPState, Array{<:Opportunity, 1}}()       # State-Action Lookup
    N = Dict{Tuple{MDPState, Opportunity}, Int32}()   # State-Action Exploration Count
    Q = Dict{Tuple{MDPState, Opportunity}, Float64}() # State-Action Reward

end

function mcts_simulate()
end

function satellite_plan_mdp_mcts(problem::PlanningProblem; sat_id::Integer=1)

    # Sort opportunities
    sort!(problem.opportunities, by = x -> x.t_start)
    
    # Set initial state
    init_opp = problem.opportunities[1]
    state = MDPState(time=init_opp.t_start, last_action=init_opp)

    plan = Opportunity[state.last_action]
    reward = 0.0

    while true
        state, action, value = mdp_mcts(problem, state)

        # Add state and action to plan
        push!(plan, state.last_action)
        reward += value

        # Break from search if terminal state
        if state.done == true
            break
        end
    end    

    return plan, reward
end