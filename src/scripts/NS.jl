cDir = @__DIR__
cd(cDir)
push!(LOAD_PATH,joinpath(cDir,"../libraries"))
#cd("C:\\Users\\rexja\\OneDrive - BYU-Idaho\\Classes\\Nelson Research Group\\MatSim.jl\\src\\libraries")
push!(LOAD_PATH,pwd())

using ase
using LennardJones
using BenchmarkTools
using YAML
using nSampling
# Load data in to test
resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")
input = YAML.load_file(joinpath(cDir,"NS.yml"))
LJ_average = LennardJones.get_LJ_averages(resultsPath)
myNS = nSampling.initialize(input["params"],LJ_average);# 10k allocations.  Need to optimize this.
nSampling.tune_step_sizes!(myNS,LJ_average)
nSampling.run_NS(myNS,LJ_average, input["output"])
T = collect(1:0.1:50000)
H,Z,Cp = nSampling.post_process(joinpath(cDir,"NS.out"),T)

plot(T[1:end],Z)
fcc = ase.fromPOSCAR(joinpath(cDir,"POSCAR.fcc"),["Ag"])
ase.eval_energy(fcc,LJ_average)