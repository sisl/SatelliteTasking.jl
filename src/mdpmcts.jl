# Exports
export MCTSState
export get_optimal_action
export mcts_rollout
export mcts_simulate
export mdp_mcts
export parallel_mdp_mcts
export satellite_plan_mdp_mcts

@with_kw mutable struct MCTSState
    T = MDPState[]                                                   # Visited States
    L = Dict{MDPState, Array{<:Opportunity, 1}}()     # State-Action Lookup
    N = Dict{Tuple{MDPState, <:Opportunity}, Int32}()   # State-Action Exploration Count
    Q = Dict{Tuple{MDPState, <:Opportunity}, Float64}() # State-Action Reward
end

"""
Get optimal action
"""
function get_optimal_action(mcts::MCTSState, state::MDPState)
    astar, vstar = :DONE, -Inf
    
    for a in mcts.L[state]
        # println("State: $(state.time) - Action: $a - Reward: $(Q[(state, a)])")
        if mcts.Q[(state, a)] > vstar
            astar, vstar = a, mcts.Q[(state, a)]
        end
    end

    return astar
end

"""
Rollout of policy
"""
function mcts_rollout(problem::PlanningProblem, mcts::MCTSState, state::MDPState, depth::Integer)

    # println("Doing Rollout from $state. Depth: $depth")

    # Exit if depth reached
    if depth == 0
        return 0.0
    end

    # Sample policy
    actions = mdp_state_actions(problem, state)

    if length(actions) == 0
        return 0.0
    end

    action_idx = ceil(Int, rand()*length(actions))
    action = actions[action_idx]

    # println("Rollout selected - Index: $action_idx - Action: $action")

    # Calculate reward for action
    reward = mdp_reward(problem, state, action)

    # Step deterministically
    next_state = mdp_step(problem, state, action)

    # Get distance action is in the future
    t_diff = abs(action.t_start - state.time)

    # Semi-markov update for rollout
    return reward + problem.solve_gamma^(t_diff)*mcts_rollout(problem, mcts, next_state, depth-1)
end

function mcts_simulate(problem::PlanningProblem, mcts::MCTSState, state::MDPState, depth::Integer)
    
    # Exit if depth reached
    if depth == 0
        return 0.0
    end

    # If state not seen, add it and rollout policy 
    if !(state in mcts.T)

        # Add state to observed states
        push!(mcts.T, state)

        # Initialize Action Lookup if not seen
        mcts.L[state] = Opportunity[]

        for action in mdp_state_actions(problem, state)

            if !(action in mcts.L[state])
                push!(mcts.L[state], action)
            end

            # Apply default rewqrd to q to encourage
            q0 = 0.0
            n0 = 0

            if typeof(action) == Sunpoint
                # pre-seed value with reward
                q0 = 1.0
                n0 = 1
            end

            # Initiallies N & Q arrays
            mcts.N[(state, action)] = n0
            mcts.Q[(state, action)] = q0
        end

        # Rollout from new state
        # println("New state, performing rollout. State: $state - Possible Actions: $(length(mcts.L[state]))")
        return mcts_rollout(problem, mcts, state, depth-1)
    end

    # println("Revisiting state, Running Simulation. State: $state")

    # Calculate N(s) by summing over all observed states 
    Ns = 0
    for a in mcts.L[state]
        Ns = Ns + mcts.N[(state, a)]
    end

    # println("Ns: $Ns")

    astar, vstar = Done(t_start=state.time), -Inf
    if Ns == 0
        # New State haven't seen this node previousy, select action randomly
        actions = mcts.L[state]
        action_idx = ceil(Int, rand()*length(actions))
        astar, vstar = actions[action_idx], Inf

        # println("Unvisited state, randomly selecting action.")
    else
        # Get optimal action given exploration bonus
        for a in mcts.L[state]
            # If already reached optimal action break iteration
            if vstar == Inf
                break
            end

            # println("MCTS Exploration: $(mcts.Q[(state, a)] + problem.mcts_c*sqrt(log(Ns)/(mcts.N[(state, a)])))")
            # print("Considering: $a - N: $(mcts.N[(state, a)])")
            if mcts.N[(state, a)] > 0 && (mcts.Q[(state, a)] + problem.mcts_c*sqrt(log(Ns)/(mcts.N[(state, a)]))) > vstar
                # print(" - Q: $(mcts.Q[(state, a)])\n")
                astar, vstar = a, mcts.Q[(state, a)]
            elseif mcts.N[(state, a)] == 0
                # Unvisited state UCB is \infty
                astar, vstar = a, Inf
                # print("\n")
            else
                # print("\n")
            end
        end
    end

    # Take optimal a
    reward = mdp_reward(problem, state, astar)
    next_state = mdp_step(problem, state, astar)

    # MCTS Simulation
    # println("Selected Action: $astar - Next State: $next_state")

    # If next state is done return reward
    if typeof(astar) == Done 
        return 0.0
    end

    # Get distance action is in the future
    t_diff = abs(astar.t_start - state.time)

    # Apply MCTS to get reward of next state
    q = reward + problem.solve_gamma^(t_diff)*mcts_simulate(problem, mcts, next_state, depth-1)
    # println("Rollout Reward: $q")

    # Update MCTS parameters 
    mcts.N[(state, astar)] = mcts.N[(state, astar)] + 1
    mcts.Q[(state, astar)] = mcts.Q[(state, astar)] + (q-mcts.Q[(state, astar)])/mcts.N[(state, astar)]

    return q
end

function mdp_mcts(problem::PlanningProblem, mcts::MCTSState, state::MDPState)
    # Run simulations from current state
    sim_outputs = []
    for i in 1:problem.mcts_sim_iterations
        # println("\nRunning new simulation")
        push!(sim_outputs, mcts_simulate(problem, mcts, state, problem.solve_depth))
    end

    # println("Simulated rewards for step: $(join([@sprintf "%4.3f" x for x in sim_outputs], ", "))")

    # Choose optimal action from rolloute and state
    action = get_optimal_action(mcts, state)

    # println("")
    # println("State: $state")
    # for a in mcts.L[state]
    #     println("Q($a): $(mcts.Q[(state, a)])")
    # end

    # Compute state action reward
    reward = mdp_reward(problem, state, action)

    if typeof(action) == Collect && action.location in state.requests
        # println("Taking duplicate collect. State: $state - Request: $(action.location.id) - MCTS: $(mcts.Q[(state,action)]) Reward: $reward")
    end

    # Step to next state with selected action
    next_state = mdp_step(problem, state, action)

    return next_state, action, reward
end

# function parallel_mdp_mcts(problem::PlanningProblem, mcts::MCTSState, state::MDPState)
#     # Run simulations from current state
#     for i in 1:problem.mcts_sim_iterations
#         q = mcts_simulate(problem, mcts, state, problem.solve_depth)
#     end

#     # Choose optimal action from rolloute and state
#     action = get_optimal_action(mcts, state)

#     # Compute state action reward
#     reward = mdp_reward(problem, state, action)

#     # Step to next state with selected action
#     state = mdp_step(problem, state, action)

#     return state, action, reward
# end

function satellite_plan_mdp_mcts(problem::PlanningProblem; parallel::Bool=false)

    # Sort opportunities
    sort!(problem.opportunities, by = x -> x.t_start)
    
    # Set initial state
    init_opp = problem.opportunities[1]
    state = MDPState(time=init_opp.t_start, last_action=init_opp)
    mcts  = MCTSState()

    states = MDPState[state]
    plan = Opportunity[state.last_action]
    action = nothing
    total_reward = 0.0
    reward  = 0.0

    while true
        if parallel == true
            # state, action, reward = parallel_mdp_mcts(problem, mcts, state)
        else
            state, action, reward = mdp_mcts(problem, mcts, state)
        end

        # println("State: $state")
        # println("Action: $action")
        # println("Actions in horizon: $(length(find_actions(problem, state, problem.solve_horizon)))")
        # println("Reward: $reward")

        # Break from search if terminal state
        if state.done == true
            # println("Returning hur. action: $action")
            break
        end

        # Add state and action to plan
        push!(states, state)
        push!(plan, action)
        total_reward += reward

        # # DEBUG: Exit Early
        # return plan, total_reward
    end    

    return states, plan, total_reward
end