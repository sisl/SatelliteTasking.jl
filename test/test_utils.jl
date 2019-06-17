let
    t_start = Epoch(2019, 1, 1, 0, 0, 0)
    mdp = MDPProblem(t_start=t_start)
    epc = Epoch(2019, 1, 1, 1, 0, 0)

    t_since = get_rel_time(mdp, epc)

    @test t_since == 3600.0
end

let
    t_start = Epoch(2019, 1, 1, 0, 0, 0)
    mdp = MDPProblem(t_start=t_start)

    t_since = 3600.0

    epc = get_abs_time(mdp, t_since)

    @test Epoch(2019, 1, 1, 1) == epc
end