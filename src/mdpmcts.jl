export satellite_plan_mdp_mcts
function satellite_plan_mdp_mcts(problem::SatPlanningProblem; enable_resources::Bool=false)

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

    # Compute initial state reward
    if typeof(state.last_action) == Collect
        r += state.last_action.location.reward
    end
    push!(plan_rewards, r)
    total_reward += r


    while !isterminal(problem, state)
        # Compute optimal action
        action = POMDPs.action(policy, state)
        
        # Advance state with next action
        next_state = POMDPs.gen(problem, state, action, Random.MersenneTwister(4)).sp
        push!(states, state)

        # Compute action reward
        r = POMDPs.reward(problem, state, action)

        # # Check next-state is safe, otherwise override with safe action
        # if next_state.power <= 0.0
        #     # Get Sunpointed action
        #     candidate_actions = problem.lt_feasible_actions[(state.last_cdo_action.id, state.last_action.id)]
        #     action = candidate_actions[findfirst(x -> typeof(x) == Sunpoint, candidate_actions)]

        #     next_state = POMDPs.gen(problem, state, action, Random.MersenneTwister(4)).sp
        #     r = POMDPs.reward(problem, state, action)
        # end

         # Alternate method for regenerating action
         if r <= 0.0 && enable_resources == true
            # Get Sunpointed action
            candidate_actions = problem.lt_feasible_actions[(state.last_cdo_action.id, state.last_action.id)]
            action = candidate_actions[findfirst(x -> typeof(x) == Sunpoint, candidate_actions)]

            next_state = POMDPs.gen(problem, state, action, Random.MersenneTwister(4)).sp

            # Compute estimated reward action
            r = POMDPs.reward(problem, state, action)
         end

        # Update State, action, reward 
        state = next_state
        push!(plan, action)
        push!(states, state)
        push!(plan_rewards, r)
        total_reward += r
    end

    return states, plan, plan_rewards, total_reward
end