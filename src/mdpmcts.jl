export satellite_plan_mdp_mcts
function satellite_plan_mdp_mcts(problem::SatPlanningProblem; parallel::Bool=false)

    # Sort opportunities
    sort!(problem.opportunities, by = x -> x.t_start)
    
    # Set initial state
    init_opp = problem.opportunities[1]
    state = initialstate(problem, Random.MersenneTwister(4))
    states = SatMDPState[state]
    plan = Opportunity[state.last_action]
    action = nothing
    total_reward = 0.0
    reward  = 0.0

    solver = MCTSSolver(n_iterations=50, depth=20, exploration_constant=5.0, enable_tree_vis=false)
    @POMDPs.requirements_info(solver, problem)
    # policy = POMDPs.solve(solver, problem);

    # state = initialstate(problem, Random.MersenneTwister(4))
    # push!(states, state)
    # reward = state.last_collect.location.reward
    # total_reward += reward

    # while state.done == false
    #     # Compute optimal action
    #     action = POMDPs.action(policy, state)
    #     push!(plan, action)
        
    #     # Advance state with next action
    #     state = POMDPs.generate_s(problem, state, action, Random.MersenneTwister(4))
    #     push!(states, state)

    #     # Update reward
    #     reward = state.last_collect.location.reward
    #     total_reward += reward
    # end

    return states, plan, total_reward
end