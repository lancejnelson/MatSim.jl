

cDir = @__DIR__
cd(cDir)
push!(LOAD_PATH,joinpath(cDir,"../libraries"))

using ase
using LennardJones
using nSampling
using YAML

# Initialize the model
resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")

# Build model with averages of the draws for fit parameters.
LJ_average = LennardJones.get_LJ_averages(resultsPath)
input = YAML.load_file(joinpath(cDir,"NS.yml"))
myNS = nSampling.initialize(input["params"],LJ_average);# Initialize the simulation...
nSampling.tune_step_sizes!(myNS,LJ_average)
nSampling.run_NS(myNS,LJ_average)


# Disregard all below

# display(ase.eval_KE(myNS.walkers[sEnergies[end-3]]))
# println(myNS.walkers[sEnergies[end-3]].coordSys)
# myNS.walkers[4].lj_vec
# ase.eval_forces(myNS.walkers[1],LJ_average)
# display(ase.get_neighbor_distances(myNS.walkers[1],4.0))
# println(myNS.walker_params.MD_time_step)
# nSampling.tune_step_sizes(myNS.walkers[sEnergies[end-2]],myNS.walker_params,LJ_average,E_max)

# ase.eval_KE(myNS.walkers[4])
# sum(0.5 * [x' * x for x in myNS.walkers[4].velocities])
# nSampling.run_NS(myNS,LJ_average)
# check =ase.get_nn_list(myNS.walkers[2],5.6)
# check = myNS.walkers[1].positions
# myvec = ["ab", "dc", "ef"]
# !("ab" in myvec)
# println(myNS.walkers[1].lVecs)
# println(ase.cell_volume(myNS.walkers[1]))
# ase.DirectToCartesian!(myNS.walkers[1])
# myNS.walker_params.MD_time_step = 0.000001
# ase.do_MD!(myNS.walkers[3],LJ_average,myNS.walker_params,E_max)
# 2 .* check #.* [2 * 5.0 ./ x for x in check] #.* 2 .* check
# [2 * 5.0 ./ x for x in check] #.* 2 .* check
# print(myNS.walkers[5].positions)
# begin_energy = ase.eval_energy(myNS.walkers[2],LJ_average)
# using Plots
# using StaticArrays
# using LinearAlgebra
# nSampling.walk_single_walker!(myNS.walkers[2],LJ_average,myNS.walker_params,ase.eval_energy(myNS.walkers[2],LJ_average))
# display(ase.eval_energy(myNS.walkers[2],LJ_average))
# ase.gradientForce(LJ_average,myNS.walkers[2],@SVector[1,4],@SVector[2,2,2])
# println(myNS.walkers[2].fitEnergy)
# display(myNS.walkers[2].lVecs)
# display(myNS.walkers[2].atomicBasis)
# ase.mapIntoCell(myNS.walkers[2])
# ase.DirectToCartesian!(myNS.walkers[2])

# using LinearAlgebra
# nHat = [[norm(x) > 0.0 ? x / norm(x) : x for x in y] for y in myNS.walkers[2].atomicBasis ]

# any([isnan(z) for x in myNS.walkers[2].atomicBasis for y in x for z in y])
# any(isnan.(myNS.walkers[2].atomicBasis[1]))
# floor(-0.5)
# rand(1:5)
# #nSampling.GMC(myNS.configs[55],500,LJ_average)
# nSampling.run_NS(myNS,LJ_average)
# rands = nSampling.get_random_displacements([2,3])

# println(rands[2])
