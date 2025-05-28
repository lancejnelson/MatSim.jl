module MatSim


currentDir = abspath(@__DIR__)
libPath = relpath(joinpath(currentDir,"libraries"))
typePath = relpath(joinpath(currentDir,"modules"))

# Add all of the public packages
#using BenchmarkTools
using LinearAlgebra
using StaticArrays
using Printf
using StatsBase
using Distributions
using DelimitedFiles
using YAML
using Plots
#using Combinatorics
# Get the types needed (MatTypes.jl in modules folder)
# nothing in here, so removing for now
# include(abspath(typePath,"MatTypes.jl"))

# Add all things that modify/handle crystal structures (CrystalUtils.jl in modules folder)
include(abspath(libPath,"crystalUtils.jl"))   # Working with Crystal structures explicitly
include(abspath(libPath,"enumeration.jl")) # Working with enumeration files to generate crystal structures or set of crystals.
include(abspath(libPath,"structuresin.jl"))  # Reading and writing dataset summary files (called structures.in)
``
# Interface with VASP (VASP.jl in modules)
include(abspath(libPath,"vaspUtils.jl"))

# Metrop Hastings Algorithms (MatML.jl in modules folder)
#include(abspath(libPath,"metrop.jl"))

# Lennard Jones model (Bundled together with MatML.jl in modules folder)
include(abspath(libPath,"Lennard-Jones/metrop.jl"))
include(abspath(libPath,"Lennard-Jones/plotting.jl"))

# Thermodynamic simulations (Not in separate module for now.)
include(abspath(libPath,"NS.jl"))


# Functions for initializing models and post processing.
include(abspath(libPath,"models.jl"))

end # module MatSim
