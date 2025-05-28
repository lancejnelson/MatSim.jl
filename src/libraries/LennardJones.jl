module LennardJones

using Distributions
#using StaticArrays
using Printf
#using Plots
using DataSets
using ase
#using LinearAlgebra
using StatsBase
#using enumeration
#using PlotsMH

# Import system libraries
using YAML
using DelimitedFiles


# The Lennard-Jones potential
struct model
    name:: String
    order:: Int64
    # params:: Array{Float64,3}
    cutoff:: Float64

    # Replace params above.
#    coeffs:: Vector{Float64}
    σ:: Vector{Float64}#T#UpperTriangular{Float64, Matrix{Float64}}
    ϵ:: Vector{Float64}#T#UpperTriangular{Float64, Matrix{Float64}}

    # Keep track of the standardization parameters so I can calculate the energy accurately.
    stdEnergy:: Float64  
    meanEnergy:: Float64
    offset:: Float64
    fitTo::String
    standardized::Bool
    system::Vector{String}
end

struct MH{D<:Distribution{Univariate,Continuous}}
    nDraws:: Int64
    nBurnIn:: Int64

    # Acceptance rates for all of the parameters
    ϵ_accept:: Vector{Float64}#Array{Float64, 2}
    σ_accept:: Vector{Float64}#Array{Float64, 2}
    std_accept:: Vector{Float64}
    
    # Prior distributions
    ϵ_Priors:: Vector{D}#Array{D,2}#Vector{Distribution}
    σ_Priors:: Vector{D}#Array{D,2}
    std_Prior:: D

    # Sigmas on proposal distributions
    ϵ_candSigs:: Vector{Float64}#Array{Float64, 2}
    σ_candSigs:: Vector{Float64}#Array{Float64, 2}
    std_candSig:: Float64

    # Initial guess
    std:: Float64


    # Proposal distribution
#    proposal:: F

    # Posterior
#    logpost:: G
end

function initialize_empty_model()
    return model("Empty",0,0.0,zeros(1),zeros(1),0.0,0.0,0.0,"None",false,[])

end

function initialize_metrop(path,model)
    metropDict = YAML.load_file(path,dicttype = Dict{String,Any})
    order = model.order
    nInteractionTypes = Int(order * (order + 1)/2)

    # Check to make sure I have all the priors needed.
    priors = metropDict["Priors"]
    if length(priors) - 1 != nInteractionTypes
        error("Number of priors specified is not consistent with the declared order")
    end


    indexDict = Dict("aa"=>[1,1], "ab"=>[1,2], "bb"=>[2,2])
    distDict = Dict("gamma"=> Gamma,"uniform"=>Uniform)
    
    # Get all of the candidate sigmas
    sigs = metropDict["candidateSigmas"]
    # Get all of the staring guesses
    start = metropDict["starting"]

   # candSig = zeros(order,order,2)
    ϵ_candSigs = zeros(nInteractionTypes)#zeros(order,order)
    σ_candSigs = zeros(nInteractionTypes)#zeros(order,order)
    for i in keys(sigs)
        if lowercase(i) != "sigma"
            index = ase.index_to_integer(indexDict[i][1],indexDict[i][2],order)
            ϵ_candSigs[index] = sigs[i]["epsilon"]
            σ_candSigs[index] = sigs[i]["sigma"]
            model.ϵ[index]= start[i]["epsilon"]
            model.σ[index]= start[i]["sigma"]
        end
    end
    std_candSig = sigs["sigma"]
    std_Guess = start["sigma"]

    # Build the array of priors
    priors = metropDict["Priors"]
    σ_Priors = Vector{Distribution{Univariate,Continuous}}(undef,nInteractionTypes)
    ϵ_Priors = Vector{Distribution{Univariate,Continuous}}(undef,nInteractionTypes)
    for i in keys(priors)
        if lowercase(i) != "sigma"
            index = ase.index_to_integer(indexDict[i][1],indexDict[i][2],order)
            σ_Priors[index]= distDict[lowercase(priors[i]["sigma"]["distribution"])](parse.(Float64,split(priors[i]["sigma"]["parameters"]))...)
            ϵ_Priors[index]= distDict[lowercase(priors[i]["epsilon"]["distribution"])](parse.(Float64,split(priors[i]["epsilon"]["parameters"]))...)
        end
    end
    std_Prior = distDict[lowercase(priors["sigma"]["distribution"])](parse.(Float64,split(priors["sigma"]["parameters"]))...)

    
    # Define the proposal distribution
    #if lowercase(metropDict["proposal"]) == "gamma"
    #end

    nDraws = metropDict["nDraws"]
    nBurnIn = metropDict["nBurnin"]
    nTotal = nBurnIn + nDraws
    ϵ_accept = zeros(nInteractionTypes)
    σ_accept = zeros(nInteractionTypes)
    std_accept = 0.0

    MH(nTotal,nBurnIn,ϵ_accept,σ_accept,[std_accept],ϵ_Priors,σ_Priors,std_Prior,ϵ_candSigs,σ_candSigs,std_candSig,std_Guess)
end



# This function takes the types of the two atoms that are interacting and returns a single integer




function initializeLJ(path::String,trainingSet;std_energy= 1.0, mean_energy = 0.0, offset = 0.0)
    modelDict = YAML.load_file(path,dicttype = Dict{String,Any})
    #println(input)
    # Initialize the LJ model     

    order = modelDict["order"]::Int64
    cutoff = modelDict["cutoff"]::Float64
    fitTo = modelDict["fitTo"]::String
    offset = modelDict["offset"]::Float64
    standardize = modelDict["standardize"]::Bool
    system = modelDict["system"]::Vector{String}
    #params = ones(order,order,2)
#    coeffs = zeros(order * (order + 1) ÷ 2)
    σ = zeros(order * (order + 1) ÷ 2)
    ϵ = zeros(order * (order + 1) ÷ 2)

    if standardize
        if fitTo == "peratom"
            mean_energy = mean([i.FP_total_energy/ase.nAtoms(i) for i in trainingSet.configs])
            std_energy = std([i.FP_total_energy/ase.nAtoms(i) for i in trainingSet.configs])
        else
            mean_energy = mean([i.FP_total_energy for i in trainingSet.configs])
            std_energy = std([i.FP_total_energy for i in trainingSet.configs])
        end

    else
        mean_energy = 0.0
        std_energy = 1.0

    end
    
    return model("LJ",order,cutoff,σ,ϵ,std_energy,mean_energy,offset,fitTo,standardize,system)
    # Pre-calculate the distances needed for LJ

#    metropDict = input["metrop"]::Dict{String,Any}
 #   LJMetrop = initializeMetrop(input["metrop"],LJ_model) 
    #return LJ_model
end





function logNormal(data::DataSets.DataSet,model,σ::Float64)::Float64
    thesum = 0.0
    for i = 1:length(data.configs)
        thesum += (data.configs[i].FP_total_energy - ase.eval_energy(data.configs[i],model))^2
    end
    thesum *= - 1/(2 * σ^2)
    
#    @time "thesum" thesum =  -data.nData/2 *log(σ^2) - 1/(2 * σ^2) * sum([(i.energyPerAtomFP - totalEnergy(i,LJ))^2   for i in data.configs])
    return -length(data.configs)/2 *log(σ^2) + thesum

end


function proposal(μ,σ)
    return Gamma(μ^2/σ^2,σ^2/μ)#proposalDict[lowercase(metrop["proposal"])]
end


function logPost(data,model,metrop,σ)::Float64
    thesum = logNormal(data,model,σ)  + logpdf(metrop.std_Prior,σ)
    for index =1:sum(1:model.order)
        thesum += logpdf(metrop.σ_Priors[index],model.σ[index])
        thesum += logpdf(metrop.ϵ_Priors[index],model.ϵ[index])
    end
    return thesum

end
function sample_σ!(metrop::MH,data::DataSets.DataSet,model,model_next,std)
    nInteractions = sum(1:model.order)
    for index=1:nInteractions
        cand = rand(proposal(model.σ[index],metrop.σ_candSigs[index]))
        if cand < 0.05
            continue
        end 
        model_next.σ[index] = cand  # Set sigma equal to the candidate draw


        # Evaluate the posterior at the candidate draw.
        numerator = logPost(data,model_next,metrop,std) + log(pdf(proposal(cand,metrop.σ_candSigs[index]),model.σ[index]))
        model_next.σ[index] = model.σ[index]  # Set sigma back to the previous draw.
    
        # Evaluate the posterior again.
        denominator =logPost(data,model_next,metrop,std)  + log(pdf(proposal(model.σ[index],metrop.σ_candSigs[index]),cand))
        r = numerator - denominator
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0.0 || ((r < 0.0) & (unif < r))  # Accept?
            model_next.σ[index] = cand
    #        metrop.params_draws[i,j,k,l] = cand   # Yes!
            metrop.σ_accept[index] += 1/metrop.nDraws
        end

    end
end



function sample_ϵ!(metrop::MH,data::DataSets.DataSet,model,model_next,std)

    for index=1:sum(1:model.order)  #Loop over all of the interaction types
        cand = rand(proposal(model.ϵ[index],metrop.ϵ_candSigs[index]))
        if cand < 0.05
            continue
        end 

        model_next.ϵ[index] = cand  # Set epsilon equal to the candidate draw

        # Evaluate the posterior at the candidate draw.
        numerator = logPost(data,model_next,metrop,std) + log(pdf(proposal(cand,metrop.ϵ_candSigs[index]),model.ϵ[index]))
        model_next.ϵ[index] = model.ϵ[index]  # Set epsilon back to the previous draw.
        
        # Evaluate the posterior again.
        denominator =logPost(data,model_next,metrop,std)  + log(pdf(proposal(model.ϵ[index],metrop.ϵ_candSigs[index]),cand))
        r = numerator - denominator
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0 || ((r < 0) & (unif < r))  # Accept?
            model_next.ϵ[index] = cand
    #        metrop.params_draws[i,j,k,l] = cand   # Yes!
            metrop.ϵ_accept[index] += 1/metrop.nDraws
        end

    end
end

function sample_std!(metrop::MH,data::DataSets.DataSet,model,prev_std)

    cand = rand(proposal(prev_std,metrop.std_candSig))
    if cand < 0.05
        return prev_std
    end 
    numerator = logPost(data,model,metrop,cand) + log(pdf(proposal(cand,metrop.std_candSig),prev_std))
    denominator =logPost(data,model,metrop,prev_std)  + log(pdf(proposal(prev_std,metrop.std_candSig),cand))
    r = numerator - denominator
    unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
    if r >= 0 || ((r < 0) & (unif < r))  # Accept?
        metrop.std_accept[1] += 1/metrop.nDraws
        return cand
        #metrop.σ_draws[i] = cand   # Yes!
    end
    return prev_std

end


function do_MH(metrop::MH,data::DataSets.DataSet,model)

    DataSets.rescaleData!(data,model.fitTo)  # Standardize the data
    #model.meanEnergy = data.meanEnergy
    #model.stdEnergy = data.stdEnergy
    #model.offset = data.offset
#

    model_next = deepcopy(model)
    intDict = Dict(1=>"a",2=>"b")
    #Write the header to the output file
    cDir = pwd()
    println("Opening file  ", joinpath(cDir,"draws.out"))
    system = "System: " * data.title * "\n"
    model_name = "Model: " * model.name * "\n"
    filename = "draws-LJ." * lstrip(data.title)
    fitTo = "fitTo: " * model.fitTo * "\n"
    io = open(joinpath(cDir,filename),"w")
    standardized = "Standardized: " * string(model.standardized) * "\n"
    mn = model.standardized ? "μ-energy: " * string(model.meanEnergy) * "\n" : "μ-energy: " * string(0.0) * "\n" 
    sn = model.standardized ? "σ-energy: " * string(model.stdEnergy) * "\n" : "σ-energy: " * string(1.0) * "\n"
    offset = model.standardized ? "offset-energy: " * string(model.offset) * "\n" : "offset-energy: " * string(0.0) * "\n"
    cutoff = "cutoff-radius: " * string(model.cutoff) * "\n"
    write(io,system)
    write(io,model_name)
    write(io,fitTo)
    write(io,standardized)
    write(io,mn)
    write(io,sn)
    write(io,offset)
    write(io,cutoff)
    acceptPosition = mark(io)  # Mark the current position
    nInteractionTypes = Int(model.order * (model.order + 1)/2)  # How many parameters do I expect to get
    nSpaces = "%" * string(15 * (2 * nInteractionTypes + 1)- 2) * "s"  # We use 15 spaces per number so let's allocate exactly the right number of spaces for the first line.
    fstring = Printf.Format(nSpaces)
    header = Printf.format(fstring, "\n")
#    header = ""
    for i =1:model.order, j = i:model.order
        header *= "         ϵ_" *intDict[i] * intDict[j] * "  "
    end

    for i =1:model.order, j = i:model.order
        header *= "         σ_" *intDict[i] * intDict[j] * "  "
    end
    header *= "         std\n"
#    header = @sprintf "ϵ_aa  ϵ_ab ϵ_bb σ_aa σ_ab σ_bb std\n"
    println(header)
    write(io,header)
    #close(io)
    std_draw = metrop.std
    for i = 1:metrop.nDraws
        println("Draw ", i)
        sample_σ!(metrop,data,model,model_next,std_draw)
        sample_ϵ!(metrop,data,model,model_next,std_draw)
        std_draw = sample_std!(metrop,data,model_next,std_draw)
        model.σ[:] .= model_next.σ[:]
        model.ϵ[:] .= model_next.ϵ[:]
        if i > metrop.nBurnIn
            writeDraw(model_next,std_draw,io)
        end
    end
    println("Acceptance rates")
    println(metrop.ϵ_accept)
    println(metrop.σ_accept)
    seek(io,acceptPosition)

    for index= 1:nInteractionTypes
        printString = @sprintf "%14.2f " metrop.ϵ_accept[index] * 100
        write(io,printString)
    end
    for index= 1:nInteractionTypes
        printString = @sprintf "%14.2f " metrop.σ_accept[index] * 100
        write(io,printString)
    end

    printString = @sprintf "%14.2f\n" metrop.std_accept[1] * 100
    write(io,printString)

    println("σ accept")
    display(metrop.σ_accept)
    println("ϵ accept")
    display(metrop.ϵ_accept)
    println("std accept")
    display(metrop.std_accept)
    close(io)
    DataSets.undoRescaling!(data,model.fitTo)  # Standardize the data

end

function writeDraw(model,std_draw,file)
    nInteractionTypes = Int(model.order * (model.order + 1)/2)  # How many parameters do I expect to get
    #io = open(file,"a")
    for index= 1:nInteractionTypes
#        index = ase.index_to_integer(i,j,model.order)
        printString = @sprintf "%15.5f" model.ϵ[index]
        write(file,printString)
    end
    for index= 1:nInteractionTypes
#        index = ase.index_to_integer(i,j,model.order)
        printString = @sprintf "%15.5f" model.σ[index]
        write(file,printString)
    end
    
    printString = @sprintf "%15.5f\n" std_draw
    write(file,printString)
end

function get_LJ_averages(filePath)
    #nInteractionTypes = Int(order * (order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    # Read the file header

    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = readHeader(filePath)
    println("cutoff")
    println(cutoff)
    #Read the draws from file
    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    #model = model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize)
    #if order != nInteractionTypes
    #    error("Order of system doesn't match with number of parameters in draw file")
    #end
    ϵ_mean = zeros(nInteractionTypes)
    σ_mean = zeros(nInteractionTypes)
    for i = 1:order, j = i:order
        for k = 1:size(data)[1]  # Loop over all the draws
            index = ase.index_to_integer(i,j,order)
    #        println(data[k,:])
            ϵ_mean[index] += convert.(Float64,data[k,(i - 1) * order + j])/nDraws
            σ_mean[index] += convert.(Float64,data[k,nInteractionTypes + (i - 1) * order + j])/nDraws
        end
    end

    return model(model_name,order, cutoff,σ_mean,ϵ_mean,sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize,system)
end


function readHeader(filePath)
    outFile = open(filePath,"r")
    system_line = split(split(readline(outFile))[2],"-")
    system = String[]
    push!(system,string(system_line[1]))
    push!(system,string(system_line[2]))
    model = split(readline(outFile))[2]
    println(system)
    println(model)
    fitTo = lowercase(split(readline(outFile))[2])
    standardize = lowercase(split(readline(outFile))[2]) == "true" ? true : false
    muEnergy = parse(Float64,split(readline(outFile))[2])
    sigmaEnergy = parse(Float64,split(readline(outFile))[2])
    offsetEnergy = parse(Float64,split(readline(outFile))[2])
    cutoff = parse(Float64,split(readline(outFile))[2])
    acceptRates = parse.(Float64,split(readline(outFile)))

    println(cutoff)
    close(outFile)
    return system,model,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates
end


end