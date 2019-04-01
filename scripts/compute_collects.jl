push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../")

# Julia Imports
using Distributions
using LinearAlgebra

# SatelliteDynamics imports
using SatelliteDynamics

# SatelliteTasking imports
using SatelliteTasking

epc0 = Epoch(2019, 1, 1, 0, 0, 0, tsys=:UTC)
epcf = Epoch(2019, 1, 2, 0, 0, 0, tsys=:UTC)

oe  = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
eci = sOSCtoCART(oe, use_degrees=true)
println(oe)
println(eci)

@time orbit = Orbit(epc0, epcf, eci, dtmax=5.0)

@time images = load_images("./data/landsat_test.json")

@time collects = find_all_collects(orbit, images, sort=true)
for opp in collects
    println(opp)
end

# oe2  = [R_EARTH + 505e3, 0, 90.0, 0, 0, 0]
# # oe2  = [R_EARTH + 500e3, 0, 90.0, 180.0, 0, 0]
# eci2 = sOSCtoCART(oe2, use_degrees=true)
# println(oe2)
# println(eci2)

pos_std = 5000.0/sqrt(3.0) # 5000 #1000 # 100
vel_std = 5.0/sqrt(3.0) # 5 # 1 # 0.1
mean     = deepcopy(eci)
covr     = diagm(0 => [pos_std^2, pos_std^2, pos_std^2, vel_std^2, vel_std^2, vel_std^2])
dist_eci = MvNormal(mean, covr)
eci2     = rand(dist_eci, 1)[:, 1]
println(eci2)
println(sCARTtoOSC(eci2, use_degrees=true))
eci_err = eci2 - eci
println("ECI Error: $eci_err")

@time orbit2 = Orbit(epc0, epcf, eci2, dtmax=5.0)
@time opportunities2 = find_all_collects(orbit2, images, sort=true)

@time opp_mean, opp_sdev = collect_stats(collects, [opportunities2])

println("Collect Mean $opp_mean")
println("Collect Std dev $opp_sdev")

# @time opp_diffs = opportunity_diff(collects, opportunities2)

# for pair in opp_diffs
#     println(pair)
# end

# error_matrix = vcat(opp_diffs...)
# println(size(vcat(opp_diffs...)))
# println(size(hcat(opp_diffs...)))
# println(error_matrix)
# println(std(error_matrix, dims=1))

# for i in 1:length(orbit.t)
#     println(orbit2.eci[:, i] - orbit.eci[:, i])
# end
