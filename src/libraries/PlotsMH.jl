module PlotsMH

using ase
using Printf
using Plots
using DataSets:DataSet
using DelimitedFiles
using LennardJones
using StatsBase
using Combinatorics
#using ase

function σ_hists(filePath)
    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = LennardJones.readHeader(filePath)
    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize)

    σ_draws = zeros(nDraws,nInteractionTypes)
    aRates = zeros(nInteractionTypes)
    for i = 1:order, j = i:order
        index = ase.index_to_integer(i,j,model.order)
        σ_draws[:,index] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
        aRates[index] = acceptRates[nInteractionTypes + i + j - 1]
    end

    println(aRates)
    
    intDict = Dict(1=>"a",2=>"b")
    x = 0:0.01:3
    σ_hists = [histogram(σ_draws[:,a],legend = false,bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("σ-%s%s\nAcceptance Rate: %5.1f %%",intDict[ase.integer_to_index(a,order)[1]],intDict[ase.integer_to_index(a,order)[2]],aRates[a]),6))) for a in 1:nInteractionTypes]

    σ = plot(σ_hists...)

    return σ
end


function ϵ_hists(filePath)

    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = LennardJones.readHeader(filePath)
    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize)

    ϵ_draws = zeros(nDraws,nInteractionTypes)
    aRates = zeros(nInteractionTypes)
    for i = 1:order, j = i:order
        index = ase.index_to_integer(i,j,model.order)
        ϵ_draws[:,index] = convert.(Float64,data[:,i + j - 1])
        aRates[index] = acceptRates[nInteractionTypes + i + j - 1]
    end
    println(aRates)
    
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    
    intDict = Dict(1=>"a",2=>"b")
    x = 0:0.01:3
    ϵ_hists = [histogram(ϵ_draws[:,a],legend = false,bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("ϵ-%s%s\nAcceptance Rate: %5.1f %%",intDict[ase.integer_to_index(a,order)[1]],intDict[ase.integer_to_index(a,order)[2]],aRates[a]),6))) for a in 1:nInteractionTypes]

    ϵ = plot(ϵ_hists...)


    return ϵ
end

function std_hist(filePath)

    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = LennardJones.readHeader(filePath)
    aRate = 0.0
    for (idx,line) in enumerate(eachline(filePath))
        if idx == 8
            aRate = parse(Float64,split(line)[end])
            break
        end
    end
 
    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize)

    std_draws = convert.(Float64,data[:,end])
    
    if convert(Int64,(size(data)[2] - 1)/2) != nInteractionTypes
        error("Order of system doesn't match with number of parameters in draw file")
    end
    

    x = 0:0.01:3
    std_hist = histogram(std_draws,legend = false,bins = 100,normalize = :pdf,annotations = ((0.5,0.95),(@sprintf("std-\nAcceptance Rate: %5.1f %%",aRate),6)))

    ϵ = plot(std_hist)
    return std_hist
end






function predPlot(filePath, dataSet::DataSet;pures::Vector{ase.atoms}= Vector{ase.atoms}(undef,2), type = "fenth")
    #Read the draws from file
    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = LennardJones.readHeader(filePath)
    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 10
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    println(system)
    model = LennardJones.model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize,system)
    ϵ_draws = zeros(nDraws,nInteractionTypes)
    σ_draws = zeros(nDraws,nInteractionTypes)
    for i = 1:order, j = i:order
        index = ase.index_to_integer(i,j,model.order)
        ϵ_draws[:,index] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,index] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end
#    if size(LJ.σ)[1] != LJ.order
#        error("Number of interaction types not matching up with specified order")
#    end

    trueVals = zeros(Float64,length(dataSet.configs))
    predictVals = zeros(Float64,length(dataSet.configs))
    predictUnc = zeros(Float64,length(dataSet.configs))
    rmsError = zeros(Float64,length(dataSet.configs))
    for j = 1:length(dataSet.configs)
        
        if fitTo == "peratom"
            trueVals[j] = dataSet.configs[j].energies[1]/dataSet.configs[j].nAtoms
        else
            trueVals[j] = dataSet.configs[j].energies[1]
        end


        overDraws = zeros(Float64,nDraws)
        println("Data point: ",j)
        for i = 1:nDraws
            model.ϵ[:] .= ϵ_draws[i,:]
            model.σ[:] .= σ_draws[i,:]
            overDraws[i] = ase.eval_energy(dataSet.configs[j],model)
        end
        predictVals[j] = mean(overDraws)
        predictUnc[j] = std(overDraws)

    end
    percentError = filter(isfinite,(predictVals .- trueVals)./trueVals)
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
    return final
end




function concentrationPlot(filePath, dataSet::DataSet;pures::Vector{ase.atoms}= Vector{ase.atoms}(undef,2), type = "fenth")
    #Read the draws from file
    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = LennardJones.readHeader(filePath)
    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize)

    ϵ_draws = zeros(nDraws,nInteractionTypes)
    σ_draws = zeros(nDraws,nInteractionTypes)
    for i = 1:order, j = i:order
        index = ase.index_to_integer(i,j,model.order)
        ϵ_draws[:,index] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,index] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end

    trueVals = zeros(Float64,length(dataSet.configs))
    predictVals = zeros(Float64,length(dataSet.configs))
    predictUnc = zeros(Float64,length(dataSet.configs))
    rmsError = zeros(Float64,length(dataSet.configs))
    for j = 1:length(dataSet.configs)
        printflag = false 
        if type == "fenth"
            trueVals[j] = ase.formationEnergy(dataSet.configs[j].energies[1]/ase.nAtoms(dataSet.configs[j]),[x.energies[1] for x in pures] ,ase.get_concentrations(dataSet.configs[j]))    
        elseif type == "peratom"
            trueVals[j] = dataSet.configs[j].energies[1]/ase.nAtoms(dataSet.configs[j])
        elseif type == "total"
            trueVals[j] = dataSet.configs[j].energies[1]
        else
            error("I don't know what kind of energies you want me to plot!")
        end

        overDraws = zeros(Float64,nDraws)
        for i = 1:nDraws
            model.ϵ[:] .= ϵ_draws[i,:]
            model.σ[:] .= σ_draws[i,:]
            overDraws[i] = ase.eval_energy(dataSet.configs[j],model)
        end
        if length(dataSet.configs[j].positions) == 1
            error("Found pures")
        end
        

        if fitTo == "peratom"
            if type == "fenth"
                if ase.get_concentrations(dataSet.configs[j])[1] == 0.0 || ase.get_concentrations(dataSet.configs[j])[2] == 0.0
                    predictVals[j] = 0.0
                    predictUnc[j] = 0.0
                else
                    predictVals[j] = mean([ase.formationEnergy(draw,[x.energies[1] for x in pures] ,ase.get_concentrations(dataSet.configs[j])) for draw in overDraws])
                    predictUnc[j] = std([ase.formationEnergy(draw,[x.energies[1] for x in pures] ,ase.get_concentrations(dataSet.configs[j])) for draw in overDraws])
                end
            elseif type == "peratom"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            elseif type == "total"
                predictVals[j] = mean(overDraws) * ase.nAtoms(dataSet.configs[j])
                predictUnc[j] = std(overDraws)* ase.nAtoms(dataSet.configs[j])
            else
                error("I don't know what kind of plot you want!")
            end
        elseif fitTo == "total"
            if type == "fenth"
                predictVals[j] = mean([ase.formationEnergy(draw/ase.nAtoms(dataSet.configs[j]),[x.energies[1] for x in pures] ,ase.get_concentrations(dataSet.configs[j])) for draw in overDraws])
                predictUnc[j] = std([ase.formationEnergy(draw/ase.nAtoms(dataSet.configs[j]),[x.energies[1] for x in pures] ,ase.get_concentrations(dataSet.configs[j])) for draw in overDraws])
            elseif type == "peratom"
                predictVals[j] = mean(overDraws)/ase.nAtoms(dataSet.configs[j])
                predictUnc[j] = std(overDraws)/ase.nAtoms(dataSet.configs[j])
            elseif type == "total"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            else
                error("I don't know what kind of plot you want!")
            end
        elseif fitTo == "fenth"
            if type == "fenth"
                predictVals[j] = mean(overDraws)
                predictUnc[j] = std(overDraws)
            elseif type == "peratom"
                predictVals[j] = mean(overDraws) + sum(dataSet.configs[j].nType/dataSet.configs[j].nAtoms .* [x.energyPerAtomFP for x in pures])
                predictUnc[j] = std(overDraws)+ sum(dataSet.configs[j].nType/dataSet.configs[j].nAtoms .* [x.energyPerAtomFP for x in pures])
            elseif type == "total"
                predictVals[j] = mean(overDraws) + sum(dataSet.configs[j].nType/dataSet.configs[j].nAtoms .* [x.energyPerAtomFP for x in pures]) * dataSet.configs[j].nAtoms
                predictUnc[j] = std(overDraws)+ sum(dataSet.configs[j].nType/dataSet.configs[j].nAtoms .* [x.energyPerAtomFP for x in pures]) * dataSet.configs[j].nAtoms
            else
                error("I don't know what kind of plot you want!")
            end

        else
            error("I don't know what kind of energies you were fitting to!")
        end
    end

    concentrations = [ase.get_concentrations(j)[1] for j in dataSet.configs]
    # Plot predicted vs true energies. Slope=1 line is a perfect fit.
    display(predictVals)
    display(trueVals)
    myp = plot(concentrations,trueVals,seriestype = :scatter,ms = 2.5,ylabel = "True Energy (eVs/atom)", xlabel = "concentration",legend=false)
    plot!(concentrations,predictVals,seriestype = :scatter,markershape = :plus,markercolor = :red,ms = 1.5,ylabel = "Predicted Energy (eVs/atom)", xlabel = "concentration",legend=false)

    savefig(myp,"cHullPlot.pdf")
    return myp
end

function tracePlots(filePath)
    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = LennardJones.readHeader(filePath)
    #Read the draws from file

    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize)

    ϵ_draws = zeros(nDraws,nInteractionTypes)
    σ_draws = zeros(nDraws,nInteractionTypes)
    for i = 1:order, j = i:order
        index = ase.index_to_integer(i,j,model.order)
        ϵ_draws[:,index] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,index] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end


    g = @layout [grid(model.order,model.order)]
    intDict = Dict(1=>"a",2=>"b")
    #,annotations = ((0.5,0.95),(@sprintf("σ-%s%s\nAcceptance Rate: %5.1f %%",intDict[a[1]],intDict[a[2]],results.ϵ_accept[a]*100),6))
    tracePlots = plot([ϵ_draws[:,a] for a in 1:nInteractionTypes]
,layout = g, size = (2000,1000),legend = false)
    plot!( [σ_draws[:,a] for a in 1:nInteractionTypes],layout = g, size = (2000,1000),legend = false)
    return tracePlots
end



function hists2d(filePath,type; ar = 1.0)
    system,model_name,fitTo,standardize,muEnergy,sigmaEnergy,offsetEnergy,cutoff,acceptRates = LennardJones.readHeader(filePath)
    #Read the draws from file

    data = readdlm(filePath,Float64;skipstart = 10)

    nDraws = countlines(filePath) - 9
    order = convert(Int64,ceil(sqrt((size(data)[2] - 1)/2)))
    nInteractionTypes = sum(1:order)
    model = LennardJones.model(model_name,order,cutoff,zeros(nInteractionTypes),zeros(nInteractionTypes),sigmaEnergy,muEnergy,offsetEnergy,fitTo,standardize)

    ϵ_draws = zeros(nDraws,nInteractionTypes)
    σ_draws = zeros(nDraws,nInteractionTypes)
    for i = 1:order, j = i:order
        index = ase.index_to_integer(i,j,model.order)
        ϵ_draws[:,index] = convert.(Float64,data[:,i + j - 1])
        σ_draws[:,index] = convert.(Float64,data[:,nInteractionTypes + i + j - 1])
    end
    #if size(model.σ)[1] != model.order
    #    error("Number of interaction types not matching up with specified order")
    #end

    intDict = Dict(1=>"a",2=>"b")
    if type == "σ-σ" || type == "ϵ-ϵ"
        combs = collect(multiset_combinations(1:nInteractionTypes,2))
        elem = [[i,j] for i = 1:model.order for j = i:model.order]
        final = [[ase.index_to_integer(elem[i[1]]...,order),ase.index_to_integer(elem[i[2]]...,order)] for i in combs]

        if type == "σ-σ"
            hist2ds = [histogram2d(σ_draws[:,x[1]...],σ_draws[:,x[2]...],xlabel = @sprintf("σ-%s%s ",intDict[ase.integer_to_index(x[1],order)[1]],intDict[ase.integer_to_index(x[1],order)[2]]),ylabel = @sprintf("σ-%s%s",intDict[ase.integer_to_index(x[2],order)[1]],intDict[ase.integer_to_index(x[2],order)[2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(σ_draws[:,x[1]...]),1.1* maximum(σ_draws[:,x[1]...])),ylim = (0.9*minimum(σ_draws[:,x[2]...]),1.1* maximum(σ_draws[:,x[2]...]))) for x in final  ]
            r = @layout [grid(length(combs),1)] 
            hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar, size = (2000,1000))
        else
            hist2ds = [histogram2d(ϵ_draws[:,x[1]...],ϵ_draws[:,x[2]...],xlabel = @sprintf("ϵ-%s%s ",intDict[ase.integer_to_index(x[1],order)[1]],intDict[ase.integer_to_index(x[1],order)[2]]),ylabel = @sprintf("ϵ-%s%s",intDict[ase.integer_to_index(x[2],order)[1]],intDict[ase.integer_to_index(x[2],order)[2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(ϵ_draws[:,x[1]...]),1.1* maximum(ϵ_draws[:,x[1]...])),ylim = (0.9*minimum(ϵ_draws[:,x[2]...]),1.1* maximum(ϵ_draws[:,x[2]...]))) for x in final  ]
            r = @layout [grid(length(combs),1)] 
            hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar, size = (2000,1000))
        end
    else
        elem = [[i,j] for i = 1:model.order for j = i:model.order]
        hist2ds = [histogram2d(σ_draws[:,ase.index_to_integer(y...,order)],ϵ_draws[:,ase.index_to_integer(x...,order)],xlabel = @sprintf("σ-%s%s ",intDict[x[1]],intDict[x[2]]),ylabel = @sprintf("ϵ-%s%s",intDict[y[1]],intDict[y[2]]),left_margin = 16Plots.mm,bottom_margin = 6Plots.mm, xlim = (0.9*minimum(σ_draws[:,ase.index_to_integer(y...,order)]),1.1* maximum(σ_draws[:,ase.index_to_integer(y...,order)])),ylim = (0.9*minimum(ϵ_draws[:,ase.index_to_integer(x...,order)]),1.1* maximum(ϵ_draws[:,ase.index_to_integer(x...,order)]))) for x in elem for y in elem ]
        r = @layout [grid(nInteractionTypes,nInteractionTypes)] 
        hist2dplots = plot(hist2ds...,layout = r,aspect_ratio = ar,colorbar=false, size = (2000,1000))

    end
    return hist2dplots
        
end

end