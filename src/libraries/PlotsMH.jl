module PlotsMH

using ase
using Printf
using Plots
using DataSets:DataSet
using DelimitedFiles
using LennardJones
using StatsBase
#using ase

function σ_hists(filePath)
    outFile = open(filePath,"r")
    system = split(readline(outFile))[1]
    standardize = lowercase(split(readline(outFile))[2]) == "true" ? true : false
    muEnergy = parse(Float64,split(readline(outFile))[2])
    sigmaEnergy = parse(Float64,split(readline(outFile))[2])
    offsetEnergy = parse(Float64,split(readline(outFile))[2])
    cutoff = parse(Float64,split(readline(outFile))[2])
    acceptRates = parse.(Float64,split(readline(outFile)))
    println(acceptRates)
    data = readdlm(filePath,Float64;skipstart = 8)

    nDraws = countlines(filePath) - 8
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(order,cutoff,zeros(order,order),zeros(order,order),sigmaEnergy,muEnergy,offsetEnergy)

    σ_draws = zeros(nDraws,order,order)
    aRates = zeros(order,order)
    for i = 1:order, j = i:order
        println(i, j,nInteractionTypes)
        println(nInteractionTypes + (i - 1) * order + j)
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
        aRates[i,j] = acceptRates[nInteractionTypes + i + j - 1]
    end
    println(aRates)
    
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    
#    if size(LJ.σ)[1] != LJ.order
#        error("Number of interaction types not matching up with specified order")
#    end

    keeps = [a for a in CartesianIndices(model.σ) if a[2] >= a[1]]
    keeps = sort(sort(keeps,by = x->x[1]),by = x->x[2])

    intDict = Dict(1=>"a",2=>"b")
    x = 0:0.01:3
    σ_hists = [histogram(σ_draws[:,a],legend = false,bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("σ-%s%s\nAcceptance Rate: %5.1f %%",intDict[a[1]],intDict[a[2]],aRates[a]),6))) for a in keeps]

    σ = plot(σ_hists...)
#    plot!(x,[pdf(results.σ_Priors[a],x) for a in keeps],lw=6,lc = :red)

    return σ
end


function ϵ_hists(filePath)

    system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff = readHeader(filePath)
    data = readdlm(filePath,Float64;skipstart = 9)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(order,cutoff,zeros(order,order),zeros(order,order),sigmaEnergy,muEnergy,offsetEnergy)

    ϵ_draws = zeros(nDraws,order,order)
    aRates = zeros(order,order)
    for i = 1:order, j = i:order
        println(i, j,nInteractionTypes)
        println(nInteractionTypes + (i - 1) * order + j)
        ϵ_draws[:,i,j] = convert.(Float64,data[:,i + j - 1])
        aRates[i,j] = acceptRates[nInteractionTypes + i + j - 1]
    end
    println(aRates)
    
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    
#    if size(LJ.ϵ)[1] != LJ.order
#        error("Number of interaction types not matching up with specified order")
#    end

    keeps = [a for a in CartesianIndices(model.ϵ) if a[2] >= a[1]]
    keeps = sort(sort(keeps,by = x->x[1]),by = x->x[2])

    intDict = Dict(1=>"a",2=>"b")
    x = 0:0.01:3
    ϵ_hists = [histogram(ϵ_draws[:,a],legend = false,bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("ϵ-%s%s\nAcceptance Rate: %5.1f %%",intDict[a[1]],intDict[a[2]],aRates[a]),6))) for a in keeps]

    ϵ = plot(ϵ_hists...)
    #plot!(x,[pdf(results.ϵ_Priors[a],x) for a in keeps],lw=6,lc = :red)


    return ϵ
end

function std_hist(filePath)

    system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff = readHeader(filePath)
    aRate = 0.0
    for (idx,line) in enumerate(eachline(filePath))
        if idx == 8
            aRate = parse(Float64,split(line)[end])
            break
        end
    end
 
    data = readdlm(filePath,Float64;skipstart = 9)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(order,cutoff,zeros(order,order),zeros(order,order),sigmaEnergy,muEnergy,offsetEnergy,fitTo)

    std_draws = convert.(Float64,data[:,end])
    
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    
#    if size(LJ.ϵ)[1] != LJ.order
#        error("Number of interaction types not matching up with specified order")
#    end

    x = 0:0.01:3
    std_hist = histogram(std_draws,legend = false,bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("std-\nAcceptance Rate: %5.1f %%",aRate),6)))

    ϵ = plot(std_hist)
#    plot!(x,pdf(results.std_Prior,x),lw =6,lc = :red)
    return std_hist
end

# " Calculates the mean Lennard Jones (LJ) parameters from a large set of data points contained within a file. Returns a LJ model object."
function LJAverages(filePath)
    #nInteractionTypes = Int(order * (order + 1)/2)  # How many parameters do I expect to get
    cDir = pwd()
    
    # Read the file header, uses key-value pairs to assign the veriables from each line
    system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff = readHeader(filePath)
    
    # Display the cutoff energy
    println("Cutoff Energy")
    println(cutoff)
    
    # Read the draws from file, accounting for whitespace
    data = readdlm(filePath,Float64;skipstart = 9)
    # Get the number of draws by subtracting the lines for the header data
    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    #if order != nInteractionTypes
    #    error("Order of system doesn't match with number of parameters in draw file")
    #end
    
    # Initialize the averages for each LJ parameter with zeros
    ϵ_mean = zeros(order,order)
    σ_mean = zeros(order,order)

    # Append the mean parameter value of each draw, normalized based on the number of draws (nDraws)
    for i = 1:order, j = i:order
        for k = 1:size(data)[1]  # Loop over all the draws
    #        println(data[k,:])
            ϵ_mean[i,j] += convert.(Float64,data[k,(i - 1) * order + j])/nDraws
            σ_mean[i,j] += convert.(Float64,data[k,nInteractionTypes + (i - 1) * order + j])/nDraws
        end
    end

    # Return the LennardJones model using the mean parameters
    return LennardJones.model(order, cutoff,σ_mean,ϵ_mean,sigmaEnergy,muEnergy,offsetEnergy,fitTo)
end

function readHeader(filePath)
    outFile = open(filePath,"r")
    system = split(readline(outFile))[1]
    fitTo = lowercase(split(readline(outFile))[2])
    standardize = lowercase(split(readline(outFile))[2]) == "true" ? true : false
    muEnergy = parse(Float64,split(readline(outFile))[2])
    sigmaEnergy = parse(Float64,split(readline(outFile))[2])
    offsetEnergy = parse(Float64,split(readline(outFile))[2])
    cutoff = parse(Float64,split(readline(outFile))[2])
    println(cutoff)
    close(outFile)
    return system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff
end

function predPlot(filePath, dataSet::DataSet;pures::Vector{ase.atoms}= Vector{ase.atoms}(undef,2), type = "fenth")
    #Read the draws from file
    system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff = readHeader(filePath)
    data = readdlm(filePath,Float64;skipstart = 9)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(order,cutoff,zeros(order,order),zeros(order,order),sigmaEnergy,muEnergy,offsetEnergy,fitTo)

    ϵ_draws = zeros(nDraws,order,order)
    σ_draws = zeros(nDraws,order,order)
    for i = 1:order, j = i:order
        println(i,j)
        println(i + j - 1)
        println(nInteractionTypes + i + j - 1)
        ϵ_draws[:,i,j] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end
#    if size(LJ.σ)[1] != LJ.order
#        error("Number of interaction types not matching up with specified order")
#    end

    trueVals = zeros(Float64,length(dataSet.crystals))
    predictVals = zeros(Float64,length(dataSet.crystals))
    predictUnc = zeros(Float64,length(dataSet.crystals))
    rmsError = zeros(Float64,length(dataSet.crystals))
    for j = 1:length(dataSet.crystals)
        if type == "fenth"
            trueVals[j] = dataSet.crystals[j].formationEnergyFP
        elseif type == "peratom"
            trueVals[j] = dataSet.crystals[j].energyPerAtomFP
        elseif type == "total"
            trueVals[j] = dataSet.crystals[j].energyPerAtomFP * dataSet.crystals[j].nAtoms
        else
            error("I don't know what kind of energies you were fitting to!")
        end

        overDraws = zeros(Float64,nDraws)
        for i = 1:nDraws
            model.ϵ[:,:] .= ϵ_draws[i,:,:]
            model.σ[:,:] .= σ_draws[i,:,:]
            overDraws[i] = LennardJones.totalEnergy(dataSet.crystals[j],model)
#            overDraws[i] = (MatSim.totalEnergy(dataSet.crystals[j],LJ) + offset) * stdEnergy + meanEnergy
        end
        if fitTo == "peratom"
            if type == "fenth"
                predictVals[j] = ase.formationEnergy(mean(overDraws),[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms) 
                predictUnc[j] = ase.formationEnergy(std(overDraws),[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms)
            elseif type == "peratom"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            elseif type == "total"
                predictVals[j] = mean(overDraws) * dataSet.crystals[j].nAtoms
                predictUnc[j] = std(overDraws)*dataSet.crystals[j].nAtoms
            else
                error("I don't know what kind of plot you want!")
            end
        elseif fitTo == "total"
            if type == "fenth"
                predictVals[j] = ase.formationEnergy(mean(overDraws)/dataSet.crystals[j]/dataSet.crystals[j].nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms)
                predictUnc[j] = ase.formationEnergy(std(overDraws)/dataSet.crystals[j]/dataSet.crystals[j].nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms)
            elseif type == "peratom"
                predictVals[j] = mean(overDraws)/dataSet.crystals[j].nAtoms
                predictUnc[j] = std(overDraws)/dataSet.crystals[j].nAtoms
            elseif type == "total"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            else
                error("I don't know what kind of plot you want!")
            end
#            predictVals[j] = ase.formationEnergy(mean(overDraws)/dataSet.crystals[j]/nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nTypes/dataSet.crystals[j].nAtoms)
#            predictUnc[j] = ase.formationEnergy(std(overDraws)/dataSet.crystals[j]/nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nTypes/dataSet.crystals[j].nAtoms)
        elseif fitTo == "fenth"
            if type == "fenth"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            elseif type == "peratom"
                predictVals[j] = mean(overDraws) + sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures])
                predictUnc[j] = std(overDraws)+ sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures])
            elseif type == "total"
                predictVals[j] = mean(overDraws) + sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures]) * dataSet.crystals[j].nAtoms
                predictUnc[j] = std(overDraws)+ sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures]) * dataSet.crystals[j].nAtoms
            else
                error("I don't know what kind of plot you want!")
            end

        else
            error("I don't know what kind of energies you were fitting to!")
        end
    end
    percentError = filter(isfinite,(predictVals .- trueVals)./trueVals)
    println(percentError)
    println(maximum(percentError))
    percentHist = histogram(percentError,bins = 100,normalize = :pdf,title = "percent error" )
    rmsError = sqrt(mean( (trueVals .- predictVals).^2 ))


    upper = maximum(trueVals)
    lower = minimum(trueVals)
    x = 1.15*lower:0.05:0.85*upper
    r = @layout [grid(2,1)]

    # Plot predicted vs true energies. Slope=1 line is a perfect fit.
    myp = plot(predictVals,trueVals,seriestype = :scatter,ms = 2.5,ylabel = "True Energy (eVs/atom)", xlabel = "Predicted Energy (eVs/atom)",legend=false)
    tag = @sprintf("RMS Error: %8.4f",rmsError)
    annotate!((0.75,0.25),tag)
    final = plot(myp,percentHist,layout = r)
    plot!(x,x,lw = 5)
    #using Printf
    return final
end




function concentrationPlot(filePath, dataSet::DataSet;pures::Vector{ase.atoms}= Vector{ase.atoms}(undef,2), type = "fenth")
    #Read the draws from file
    system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff = readHeader(filePath)
    data = readdlm(filePath,Float64;skipstart = 9)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(order,cutoff,zeros(order,order),zeros(order,order),sigmaEnergy,muEnergy,offsetEnergy,fitTo)

    ϵ_draws = zeros(nDraws,order,order)
    σ_draws = zeros(nDraws,order,order)
    for i = 1:order, j = i:order
        ϵ_draws[:,i,j] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end
#    if size(LJ.σ)[1] != LJ.order
#        error("Number of interaction types not matching up with specified order")
#    end

    trueVals = zeros(Float64,length(dataSet.crystals))
    predictVals = zeros(Float64,length(dataSet.crystals))
    predictUnc = zeros(Float64,length(dataSet.crystals))
    rmsError = zeros(Float64,length(dataSet.crystals))
    for j = 1:length(dataSet.crystals) 
        if type == "fenth"
            trueVals[j] = dataSet.crystals[j].formationEnergyFP
        elseif type == "peratom"
            trueVals[j] = dataSet.crystals[j].energyPerAtomFP
        elseif type == "total"
            trueVals[j] = dataSet.crystals[j].energyPerAtomFP * dataSet.crystals[j].nAtoms
        else
            error("I don't know what kind of energies you were fitting to!")
        end

        overDraws = zeros(Float64,nDraws)
        for i = 1:nDraws
            model.ϵ[:,:] .= ϵ_draws[i,:,:]
            model.σ[:,:] .= σ_draws[i,:,:]
            overDraws[i] = LennardJones.totalEnergy(dataSet.crystals[j],model)
#            overDraws[i] = (MatSim.totalEnergy(dataSet.crystals[j],LJ) + offset) * stdEnergy + meanEnergy
        end
        if fitTo == "peratom"
            if type == "fenth"
                predictVals[j] = ase.formationEnergy(mean(overDraws),[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms) 
                predictUnc[j] = ase.formationEnergy(std(overDraws),[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms)
            elseif type == "peratom"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            elseif type == "total"
                predictVals[j] = mean(overDraws) * dataSet.crystals[j].nAtoms
                predictUnc[j] = std(overDraws)*dataSet.crystals[j].nAtoms
            else
                error("I don't know what kind of plot you want!")
            end
        elseif fitTo == "total"
            if type == "fenth"
                predictVals[j] = ase.formationEnergy(mean(overDraws)/dataSet.crystals[j]/dataSet.crystals[j].nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms)
                predictUnc[j] = ase.formationEnergy(std(overDraws)/dataSet.crystals[j]/dataSet.crystals[j].nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms)
            elseif type == "peratom"
                predictVals[j] = mean(overDraws)/dataSet.crystals[j].nAtoms
                predictUnc[j] = std(overDraws)/dataSet.crystals[j].nAtoms
            elseif type == "total"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            else
                error("I don't know what kind of plot you want!")
            end
#            predictVals[j] = ase.formationEnergy(mean(overDraws)/dataSet.crystals[j]/nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nTypes/dataSet.crystals[j].nAtoms)
#            predictUnc[j] = ase.formationEnergy(std(overDraws)/dataSet.crystals[j]/nAtoms,[x.energyPerAtomFP for x in pures] ,dataSet.crystals[j].nTypes/dataSet.crystals[j].nAtoms)
        elseif fitTo == "fenth"
            if type == "fenth"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            elseif type == "peratom"
                predictVals[j] = mean(overDraws) + sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures])
                predictUnc[j] = std(overDraws)+ sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures])
            elseif type == "total"
                predictVals[j] = mean(overDraws) + sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures]) * dataSet.crystals[j].nAtoms
                predictUnc[j] = std(overDraws)+ sum(dataSet.crystals[j].nType/dataSet.crystals[j].nAtoms .* [x.energyPerAtomFP for x in pures]) * dataSet.crystals[j].nAtoms
            else
                error("I don't know what kind of plot you want!")
            end

        else
            error("I don't know what kind of energies you were fitting to!")
        end
    end
#    percentError = (predictVals .- trueVals)./trueVals
#    percentHist = histogram(percentError,bins = 100,normalize = :pdf,title = "percent error")
#    rmsError = sqrt(mean( (trueVals .- predictVals).^2 ))


 #   upper = maximum(trueVals)
 #   lower = minimum(trueVals)
 #   x = 1.15*lower:0.05:0.85*upper
 #   r = @layout [grid(2,1)]
    concentrations = [(j.nType/j.nAtoms)[1] for j in dataSet.crystals]
    # Plot predicted vs true energies. Slope=1 line is a perfect fit.
    display(predictVals)
    display(trueVals)
    myp = plot(concentrations,trueVals,seriestype = :scatter,ms = 2.5,ylabel = "True Energy (eVs/atom)", xlabel = "concentration",legend=false)
    plot!(concentrations,predictVals,seriestype = :scatter,markershape = :plus,markercolor = :red,ms = 1.5,ylabel = "Predicted Energy (eVs/atom)", xlabel = "concentration",legend=false)
 #   tag = @sprintf("RMS Error: %8.4f",rmsError)
 #   annotate!((0.75,0.25),tag)
 #   final = plot(myp,percentHist,layout = r)
 #   plot!(x,x,lw = 5)
    #using Printf
    savefig(myp,"cHullPlot.pdf")
    return myp
end

function tracePlots(filePath)
    system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff = readHeader(filePath)
    #Read the draws from file

    data = readdlm(filePath,Float64;skipstart = 9)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(order,cutoff,zeros(order,order),zeros(order,order),sigmaEnergy,muEnergy,offsetEnergy,fitTo)

    ϵ_draws = zeros(nDraws,order,order)
    σ_draws = zeros(nDraws,order,order)
    for i = 1:order, j = i:order
        ϵ_draws[:,i,j] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end
    if size(model.σ)[1] != model.order
        error("Number of interaction types not matching up with specified order")
    end
    # We don't use the lower left triangle of the matrix of parameters, so let's get the indices right now and sort them so the plots
    # are arranged correctly.
    keeps = [a for a in CartesianIndices(model.σ) if a[2] >= a[1]]
    keeps = sort(sort(keeps,by = x->x[1]),by = x->x[2])

    g = @layout [grid(model.order,model.order)]
    intDict = Dict(1=>"a",2=>"b")
    #,annotations = ((0.5,0.95),(@sprintf("σ-%s%s\nAcceptance Rate: %5.1f %%",intDict[a[1]],intDict[a[2]],results.ϵ_accept[a]*100),6))
    tracePlots = plot([ϵ_draws[:,a] for a in keeps]
,layout = g, size = (2000,1000),legend = false)
    plot!( [σ_draws[:,a] for a in keeps],layout = g, size = (2000,1000),legend = false)
    return tracePlots
end



function hists2d(filePath,type; ar = 1.0)
    system,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff = readHeader(filePath)
    #Read the draws from file

    data = readdlm(filePath,Float64;skipstart = 9)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(order,cutoff,zeros(order,order),zeros(order,order),sigmaEnergy,muEnergy,offsetEnergy,fitTo)

    ϵ_draws = zeros(nDraws,order,order)
    σ_draws = zeros(nDraws,order,order)
    for i = 1:order, j = i:order
        ϵ_draws[:,i,j] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,i,j] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end
    if size(model.σ)[1] != model.order
        error("Number of interaction types not matching up with specified order")
    end

    intDict = Dict(1=>"a",2=>"b")
    if type == "σ-σ" || type == "ϵ-ϵ"
        combs = collect(multiset_combinations(1:nInteractionTypes,2))
        elem = [[i,j] for i = 1:model.order for j = i:model.order]
        final = [[elem[i[1]],elem[i[2]]] for i in combs]

        if type == "σ-σ"
            hist2ds = [histogram2d(σ_draws[:,x[1]...],σ_draws[:,x[2]...],xlabel = @sprintf("σ-%s%s ",intDict[x[1][1]],intDict[x[1][2]]),ylabel = @sprintf("σ-%s%s",intDict[x[2][1]],intDict[x[2][2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(σ_draws[:,x[1]...]),1.1* maximum(σ_draws[:,x[1]...])),ylim = (0.9*minimum(σ_draws[:,x[2]...]),1.1* maximum(σ_draws[:,x[2]...]))) for x in final  ]
            r = @layout [grid(length(combs),1)] 
            hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar, size = (2000,1000))
        else
            hist2ds = [histogram2d(ϵ_draws[:,x[1]...],ϵ_draws[:,x[2]...],xlabel = @sprintf("ϵ-%s%s ",intDict[x[1][1]],intDict[x[1][2]]),ylabel = @sprintf("ϵ-%s%s",intDict[x[2][1]],intDict[x[2][2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(ϵ_draws[:,x[1]...]),1.1* maximum(ϵ_draws[:,x[1]...])),ylim = (0.9*minimum(ϵ_draws[:,x[2]...]),1.1* maximum(ϵ_draws[:,x[2]...]))) for x in final  ]
            r = @layout [grid(length(combs),1)] 
            hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar, size = (2000,1000))
        end
    else
        elem = [[i,j] for i = 1:model.order for j = i:model.order]
        hist2ds = [histogram2d(σ_draws[:,y...],ϵ_draws[:,x...],xlabel = @sprintf("σ-%s%s ",intDict[x[1]],intDict[x[2]]),ylabel = @sprintf("ϵ-%s%s",intDict[y[1]],intDict[y[2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(σ_draws[:,y...]),1.1* maximum(σ_draws[:,y...])),ylim = (0.9*minimum(ϵ_draws[:,x...]),1.1* maximum(ϵ_draws[:,x...]))) for x in elem for y in elem ]
        r = @layout [grid(nInteractionTypes,nInteractionTypes)] 
        hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar,colorbar=false, size = (2000,1000))

    end
    return hist2dplots
        
end

end