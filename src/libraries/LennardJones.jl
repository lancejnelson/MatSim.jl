module LennardJones

# Import local libraries
#using Distributions
#using StaticArrays
#using Printf
#using Plots
#using DataSets
#using ase
#using LinearAlgebra
#using enumeration
#using PlotsMH

# Import system libraries
using YAML


# The Lennard-Jones potential
struct model{T}
    order:: Int64
    # params:: Array{Float64,3}
    cutoff:: Float64

    # Replace params above.
    # sigma is the atomic spacing parameter in nm
    σ:: T#UpperTriangular{Float64, Matrix{Float64}}
    # epsilon when divided by k_B, Boltzman's constant, is a temperature in Kelvin
    ϵ:: T#UpperTriangular{Float64, Matrix{Float64}}

    # Keep track of the standardization parameters so I can calculate the energy accurately.
    stdEnergy:: Float64  
    meanEnergy:: Float64
    offset:: Float64
    fitTo::String

end


function initializeLJ(path::String)
    modelDict = YAML.load_file(path,dicttype = Dict{String,Any})
    #println(input)
    # Initialize the LJ model     

    order = modelDict["order"]::Int64
    cutoff = modelDict["cutoff"]::Float64
    fitTo = modelDict["fitTo"]::String
    #params = ones(order,order,2)
    σ = zeros(order,order)
    ϵ = zeros(order,order)
    
    return model(order,cutoff,σ,ϵ,1.0,0.0,0.0,fitTo)
    # Pre-calculate the distances needed for LJ

#    metropDict = input["metrop"]::Dict{String,Any}
 #   LJMetrop = initializeMetrop(input["metrop"],LJ_model) 
    return LJ_model
end










end