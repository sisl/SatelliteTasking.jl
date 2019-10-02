export satellite_plan_mdp_mcts
function satellite_plan_mdp_mcts(problem::SatPlanningProblem; parallel::Bool=false)

    # Sort opportunities
    sort!(problem.opportunities, by = x -> x.t_start)
    
    # Set initial state
    init_opp = problem.opportunities[1]
    state = initialstate(problem, Random.MersenneTwister(4))
    states = SatMDPState[state]
    plan = Opportunity[state.last_action]
    plan_rewards = Real[]
    action = nothing
    total_reward = 0.0
    r  = 0.0

    solver = MCTSSolver(n_iterations=problem.mcts_n_iterations, depth=problem.solve_depth, exploration_constant=problem.mcts_exploration_constant)
    # @POMDPs.requirements_info(solver, problem)
    policy = POMDPs.solve(solver, problem);

    state = initialstate(problem, Random.MersenneTwister(4))
    push!(states, state)
    # r = state.last_collect.location.r
    # r = reward(problem, state, state.last_action)
    if typeof(state.last_action) == Collect
        r = state.last_action.location.reward
    else
        r = 0
    end
    push!(plan_rewards, r)
    total_reward += r

    # println("State: $state\n")
    # println("Reward: $r\n")

    while !isterminal(problem, state)
        # Compute optimal action
        action = POMDPs.action(policy, state)
        push!(plan, action)
        # println("Action: $action\n")
        
        # Advance state with next action
        state = POMDPs.gen(problem, state, action, Random.MersenneTwister(4)).sp
        push!(states, state)
        # println("State: $state\n")

        # Update r
        # r = state.last_collect.location.r
        if typeof(state.last_action) == Collect
            r = state.last_action.location.reward
        else
            r = 0
        end
        push!(plan_rewards, r)
        # println("Request Ids: $(state.request_ids)")
        # println("Reward: $r\n")
        total_reward += r
    end

    return states, plan, plan_rewards, total_reward
end