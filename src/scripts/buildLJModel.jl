cd("/Users/legoses/OneDrive - BYU-Idaho/codes/MatSim/src/libraries/")
push!(LOAD_PATH,pwd())


using SOAP
using ase
using YAML
using enumeration
using LennardJones
using PlotsMH
using nSampling
using StillingerWeber
using DataSets
using vaspUtils

cDir = @__DIR__
cd(cDir)
# Get the training and holdout sets
trainingSet,holdoutSet = DataSets.getTraining_Holdout_Sets(joinpath(cDir,"structures.AgPt"),species = ["Pt","Ag"]);

# Initialize the model
LJ = LennardJones.initializeLJ(joinpath(cDir,"LJ.yml"),trainingSet);
# Initialize the MH setting
MH = LennardJones.initialize_metrop(joinpath(cDir,"metrop.yml"),LJ)
# Pre-calculate the distances needed for LJ
DataSets.precalc_LJ_distances!(trainingSet,LJ.cutoff)
DataSets.precalc_LJ_distances!(holdoutSet,LJ.cutoff)
# Get samples
LennardJones.do_MH(MH,trainingSet,LJ)

resultsPath = joinpath(cDir,"draws-LJ.Pt-Ag")

# Plot the results
pureCrystals = DataSets.findPureEnergies(joinpath(cDir,"structures.AgPt"))
PlotsMH.concentrationPlot(resultsPath,holdoutSet,pures=pureCrystals,type = "fenth")
PlotsMH.predPlot(resultsPath,holdoutSet,pures=pureCrystals)
PlotsMH.std_hist(resultsPath)
hists = PlotsMH.std_hist(resultsPath)
PlotsMH.tracePlots(resultsPath)
PlotsMH.σ_hists(resultsPath)
PlotsMH.ϵ_hists(resultsPath)

PlotsMH.hists2d(resultsPath,"ϵ-σ";ar = :none)
# Build model with averages of the draws for fit parameters.
LJ_average = PlotsMH.LJAverages(resultsPath)
species = split(trainingSet.title,"-")
path = joinpath(cDir,"struct_enum.out")
DataSets.gss(path,LJ_average,species,readend = 10800)  # Predict energies of everything in struct_enum.out