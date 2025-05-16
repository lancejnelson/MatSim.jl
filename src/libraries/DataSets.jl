module DataSets

using ase
using StatsBase
using Printf
using YAML
#using metrop
#using LazySets
using QHull
#using LennardJones
using enumeration

struct DataSet
    title:: String
    configs::Vector{ase.atoms}
#    stdEnergy:: Float64
#    meanEnergy::Float64
#    offset::Float64
#    nData::Int
#    standardized::Bool
#    fitTo::String
end


function initialize(file;species = ["N-A", "N-A"],offset = 0.0)
    # If you give a structures.in file directly, read directly from that. 
    if occursin("structures",file)
        println("Reading from" * file)
        dset = fromStructuresIn(file,species,overwriteLatPar = false,offset = offset)
        # Add other test cases as needed
    # Otherwise, parse the YAML file and read from the appropriate file.
    elseif occursin("yml",file)
        input = YAML.load_file(file,dicttype = Dict{String,Any})

        
        # Get the dataset
        species = String[x for x in input["dataset"]["species"]]
        dataFile = input["dataset"]["file"]
        println("Reading from" * dataFile)
        dset = fromStructuresIn(dataFile,species,overwriteLatPar = false,offset = offset)
    
    end
    return dset
end









function gss(file,model,species;readend = 100)
#    filePath = joinpath(dirname(file), "structures.AgPt")
#    pures = findPureEnergies(filePath)
    pureatoms = ase.fccPures(species)
    pures = []
    for pure in pureatoms
        push!(pures,ase.eval_energy(pure,model))
    end
    enum=enumeration.read_header(file)

    cDir = dirname(file)
    io = open(joinpath(cDir,"gss.out"),"w")
#    for (idx,line) in enumerate(eachline(file))
    for idx in 1:readend
        #println(idx)
 #       if idx > readend
 #           break
 #       end
       # if idx < 16 
       #     continue
       # end



#        hnfN = parse(Int,split(line)[2])
#        hnf_degen = parse(Int,split(line)[3])
#        label_degen = parse(Int,split(line)[4])
#        total_degen = parse(Int,split(line)[5])
#        sizeN = parse(Int,split(line)[6])
#        n = parse(Int,split(line)[7])
#        pgOps = parse(Int,split(line)[8])
#        SNF = Diagonal(parse.(Int,split(line)[9:11]))
#        a = parse(Int,split(line)[12])
#        b = parse(Int,split(line)[13])
#        c = parse(Int,split(line)[14])
#        d = parse(Int,split(line)[15])
#        e = parse(Int,split(line)[16])
#        f = parse(Int,split(line)[17])
#        HNF = LowerTriangular([a 0 0
#                               b c 0 
#                               d e f])
#
#        l = parse.(Int,split(line)[18:26])
#        lTransform = hcat([l[i:i+2] for i=1:3:7]...)'
#        println(split(line))
#        labeling = split(line)[27]
#        arrows = try split(line)[28] catch y repeat("0",length(labeling)) end
#        strN = idx - 15
#        eStruct =  EnumStruct(strN,hnfN,hnf_degen,label_degen,total_degen,sizeN,n,pgOps,SNF,HNF,lTransform,labeling,arrows)
#        println(idx)
        atoms = ase.fromEnum(file,idx,String.(species))
#        atoms = ase.fromEnum(file,idx,["Na","Na"])
        #display(atoms)

#        atoms = ase.atoms(enum,eStruct,["Ag","Pt"],mink=true)
        atoms.model_energy = ase.eval_energy(atoms,model)
        conc = ase.get_concentrations(atoms)
        formationEnergyModel = ase.formationEnergy(atoms.model_energy,pures,conc)
        printString = @sprintf "%5d  %8.4f %8.4f %8.4f %8.4f\n" idx conc[1] conc[2] atoms.model_energy formationEnergyModel
        write(io,printString)
    end
    close(io)
end







#DataSet(dSet) = fromases(dSet,standardize)

function getConvexHull(file::String)

    nLines = countlines(file)  # Count the number of lines in the file.
    data = zeros(nLines,2)
    foundPureA = false
    foundPureB = false
    for (idx,line) in enumerate(eachline(file))  # Loop over each line in the file.
        concs = parse.(Float64,split(line)[2:3])  # Get the concentrations
        if isapprox(concs,[0.0,1.0])
            foundPureA = true
        elseif isapprox(concs,[1.0,0.0])
            foundPureB = true
        end
        data[idx,:] = parse.(Float64,split(line)[2:3:5])  # Get the second (concentration) and fifth (formation energy) entries in the file.
    #    push!(data,parse.(Float64,split(line)[2:3:5]))  # Get the second (concentration) and fifth (formation energy) entries in the file.
    end

    if !foundPureA
        data = vcat(data,[0.0 0.0])  # Add the pure A point if it is not present.
    end
    
    if !foundPureA
        data = vcat(data,[1.0 0.0])  # Add the pure A point if it is not present.
    end
    hull = QHull.chull(data)  # Find the convex hull of the data points.
    positiveRemoved = zeros(length(hull.vertices),2)  #Initialize an empty vector to hold the points that are below the zero formation energy line.
    counter = 1
    for (idx,value) in enumerate(hull.vertices) # Loop over the hull and only keep the points that are below the zero formation energy line.
        println(value)
        println(hull.points[value])
        if hull.points[value,2] <= 0.0 
            positiveRemoved[counter,:] = hull.points[value,:]  # Store the points that are below the zero formation energy line.
            counter += 1
        end
    end
    perms = sortperm(positiveRemoved[:,1])  # Sort the points by concentration.
    positiveRemoved = positiveRemoved[perms,:]  # Reorder the points by concentration.
    return unique(positiveRemoved,dims =1)  # Since I started with a matrix of zeros and then removed points, there may be unnecessary zeros remaining in the array. Remove them.
end

function undoRescaling!(dataSet::DataSet,fitTo)
    #meanEnergy = mean([i.FP_total_energy for i in dataSet.configs])
    #stdEnergy = std([i.FP_total_energy for i in dataSet.configs])
    #offset = 3
    for i in dataSet.configs
        if fitTo == "peratom"
            i.FP_total_energy = i.FP_total_energy * ase.nAtoms(i)
#        else
#            i.FP_total_energy = (i.FP_total_energy + dataSet.offset)*dataSet.stdEnergy+dataSet.meanEnergy
        end
    end
    #dataSet.offset = offset
#    dataSet.standardized = false
    return nothing #stdEnergy, meanEnergy
end
function standardizeData!(dataSet::DataSet)
    #meanEnergy = mean([i.FP_total_energy for i in dataSet.configs])
    #stdEnergy = std([i.FP_total_energy for i in dataSet.configs])
    #offset = 3
    for i in dataSet.configs
        i.FP_total_energy = (i.FP_total_energy - dataSet.meanEnergy)/dataSet.stdEnergy-dataSet.offset
    end
    #dataSet.offset = offset
#    dataSet.standardized = true
    return nothing #stdEnergy, meanEnergy
end
#function set_training_set(dSet,standardize,fitTo)
#    if standardize
#        println("Standardizing Data")
#        if lowercase(fitTo) == "peratom"
#            println("Using per atom energies")
#            for i in dSet.configs
#                i.fitEnergy = (i.FP_total_energy - dSet.meanEnergy)/dSet.stdEnergy-dSet.offset
#            end
#        elseif lowercase(fitTo) == "total"
#            println("Using total energies")
#            for i in dSet.configs
#                i.fitEnergy = (i.FP_total_energy * ase.nAtoms(i) - dSet.meanEnergy)/dSet.stdEnergy-dSet.offset
#            end
#        elseif lowercase(fitTo) == "fenth"
#            println("Using formation energies")
#            for i in dSet.configs
#                i.fitEnergy = (i.formationEnergyFP - dSet.meanEnergy)/dSet.stdEnergy-dSet.offset
#            end
#        end
#        return DataSet(dSet.title,dSet.configs,dSet.stdEnergy,dSet.meanEnergy,dSet.offset,dSet.nData,true,lowercase(fitTo))
#    else
#        println("Not standardizing data")
#        if lowercase(fitTo) == "peratom"
#            println("Using per atom energies")
#            for i in dSet.configs
#                i.fitEnergy = i.FP_total_energy
#            end
#        elseif lowercase(fitTo) == "total"
#            println("Using total energies")
#            for i in dSet.configs
#                i.fitEnergy = i.FP_total_energy * ase.nAtoms(i)
#            end
#        elseif lowercase(fitTo) == "fenth"
#            println("Using formation energies")
#            for i in dSet.configs
#                i.fitEnergy = i.formationEnergyFP
#            end
#        end
#        return DataSet(dSet.title,dSet.configs,dSet.stdEnergy,dSet.meanEnergy,dSet.offset,dSet.nData,false,lowercase(fitTo))
#
#    end
#
#end

function getTraining_Holdout_Sets(file::String;nTraining=100,species = ["N-A", "N-A"])
    # If you give a structures.in file directly, read directly from that. 
    if occursin("structures",file)
        println("Reading from" * file)
        dset = fromStructuresIn(file,species,overwriteLatPar = false)
        # Add other test cases as needed
    # Otherwise, parse the YAML file and read from the appropriate file.
    elseif occursin("yml",file)
        input = YAML.load_file(file,dicttype = Dict{String,Any})

        
        # Get the dataset
        species = String[x for x in input["dataset"]["species"]]
        dataFile = input["dataset"]["file"]
        offset = Float64(input["dataset"]["offset"])
        standardize = Bool(input["dataset"]["standardize"])
        fitTo = String(input["dataset"]["fitTo"])
        nTraining = Int(input["dataset"]["nTraining"])

        println("Reading from" * dataFile)
        dset = fromStructuresIn(dataFile,species,overwriteLatPar = false)

    end
    trainingSet,holdoutSet = getTraining_Holdout_Sets(dset,nTraining)

    return trainingSet,holdoutSet
end


function precalc_LJ_distances!(dset::DataSet,cutoff)
    for i in dset.configs
        ase.precalc_LJ_distances!(i,cutoff)
    end

end

function nConfigs(dset::DataSet)
    return length(dset.configs)
end
function getTraining_Holdout_Sets(dset::DataSet,nStructures)

    training = StatsBase.sample(1:length(dset.configs),nStructures,replace = false)
    holdout = setdiff(1:length(dset.configs),training)

#    if fitTo == "peratom"
#        meanTraining = mean([i.FP_total_energy/ase.nAtoms(i) for i in dset.configs[training]])
#        stdTraining = std([i.FP_total_energy/ase.nAtoms(i) for i in dset.configs[training]])
#        meanHoldout = mean([i.FP_total_energy/ase.nAtoms(i) for i in dset.configs[holdout]])
#        stdHoldout = std([i.FP_total_energy/ase.nAtoms(i) for i in dset.configs[holdout]])
#    else
#        meanTraining = mean([i.FP_total_energy for i in dset.configs[training]])
#        stdTraining = std([i.FP_total_energy for i in dset.configs[training]])
#        meanHoldout = mean([i.FP_total_energy for i in dset.configs[holdout]])
#        stdHoldout = std([i.FP_total_energy for i in dset.configs[holdout]])
#    end

    trainingSet = DataSet(dset.title,dset.configs[training])
    holdoutSet = DataSet(dset.title,dset.configs[holdout])

    return trainingSet, holdoutSet#set_training_set(trainingSet,standardize,fitTo), holdoutSet
    
end

function findPureEnergies(filePath)
    # Get the order of the system from the first POSCAR file.
    order = 0
    for (idx,line) in enumerate(eachline(filePath))
        if idx == 8
            order = length(split(line))
            break
        end
        if idx > 100
            error("We shouldn't have to read this far into the structures.in file to find the order")

        end
    end
    pures = Vector{ase.atoms}(undef,order)  # Initialize a vector to hold the pure substances.
    cLines = Vector{String}()
    counter = 0
    energyType = "peratom"
    for (idx,line) in enumerate(eachline(filePath))
        if idx == 1
            energyType = split(line)[1]
            continue
        end
        if idx < 3
            continue
        end
        if occursin("#--",line)
            thisCrystal = ase.fromPOSCAR(cLines,["N-a", "N-a"])
#            energy = parse(Float64,cLines[end])
            if lowercase(energyType) == "peratom"
                thisCrystal.FP_total_energy = parse(Float64,cLines[end]) * ase.nAtoms(thisCrystal)
            elseif lowercase(energyType) == "total"
                thisCrystal.FP_total_energy = parse(Float64,cLines[end])
#            elseif lowercase(energyType) == "fenth"
#                thisCrystal.formationEnergyFP = parse(Float64,cLines[end])
            else 
                error("Can't recognize the first line of structures.in")
            end
            if count(!iszero,ase.nType(thisCrystal)) == 1  # If there is only one nonzero entry in the nType vector, it's a pure configuration.
                pures[argmax(ase.nType(thisCrystal))] = thisCrystal
            end
            cLines = Vector{String}()
            counter = 0
        else
            push!(cLines,line)
            counter += 1
        end
    end
    return pures
end
function fromStructuresIn(filePath,species::Vector{String};overwriteLatPar = false,offset = 0.0)
    # Get the energies of the pure substances first so we can calculate formation energies
#    pures = findPureEnergies(filePath)
#    foundPures = true
#    if all(pures .== 0)
#        foundPures = false
#        println(" Couldn't find pure energies. Formation Energies will not be accurate")
#    end

    # The species list should always be in reverse alphabetical order. (All VASP calculations are performed under that convention)
    sort!(species,rev = true)

    data = Vector{ase.atoms}()
    cLines = Vector{String}()
    title = join(species,"-")
    counter = 0
    energyType = "peratom"
    for (idx,line) in enumerate(eachline(filePath))
        # The first line should tell me what kind of energies are present.
        if idx == 1
            energyType = split(line)[1]
            continue
        end
        if idx < 3
            continue
        end
        if occursin("#--",line)
            thisCrystal = ase.fromPOSCAR(cLines,species,overwriteLatPar = overwriteLatPar)
            if lowercase(energyType) == "peratom"
                thisCrystal.FP_total_energy = parse(Float64,cLines[end]) * ase.nAtoms(thisCrystal)
#                if foundPures
#                    println("calculating formation energy!")
#                    println(thisCrystal.formationEnergyFP)
#                    concentrations = ase.nType(thisCrystal) /ase.nAtoms(thisCrystal)
#                    thisCrystal.formationEnergyFP = ase.formationEnergy(thisCrystal.FP_total_energy,[x.FP_total_energy/ase.nAtoms(x) for x in pures],concentrations)
#                    MatSim.formationEnergy!(thisCrystal,pures)
#                else
#                    thisCrystal.formationEnergyFP = NaN
 #               end
            elseif lowercase(energyType) == "total"
                thisCrystal.FP_total_energy = parse(Float64,cLines[end])
#                if foundPures
#                    concentrations = ase.nType(thisCrystal) /ase.nAtoms(thisCrystal)
#                    thisCrystal.formationEnergyFP = ase.formationEnergy(thisCrystal.FP_total_energy,[x.FP_total_energy/ase.nAtoms(x) for x in pures],concentrations)
##                    MatSim.formationEnergy!(thisCrystal,pures)
#                else
#                    thisCrystal.formationEnergyFP = NaN
#                end
#            elseif lowercase(energyType) == "fenth"
#                thisCrystal.formationEnergyFP = parse(Float64,cLines[end])
#                if foundPures
#                    ase.totalEnergyFromFormationEnergy!(thisCrystal,[x.FP_total_energy/ase.nAtoms(x) for x in pures])
#                else
#                    thisCrystal.FP_total_energy = NaN
#                end
            else 
                error("Can't recognize the first line of structures.in")
            end
            if !isnan(thisCrystal.FP_total_energy)
                push!(data,thisCrystal)
            end
            cLines = Vector{String}()
            counter = 0
        else
            push!(cLines,line)
            counter += 1
        end
#        if occursin("#--",line)
#            nAtoms = sum(parse(Int64,x) for x in split(pos[idx + 6]))
#            startpoint = idx + 1
#            theend = idx + 7 + nAtoms
#            thisCrystal = ase(pos[startpoint:theend],species,overwriteLatPar = overwriteLatPar, FP_total_energy = parse(Float64,pos[theend + 2 ]))
#            if !isnan(thisCrystal.FP_total_energy)
#                push!(data,thisCrystal)
#            end
#        end
    end
#    meanEnergy = mean([i.FP_total_energy for i in data])
#    stdEnergy = std([i.FP_total_energy for i in data])    
#    println("Mean Energy: ",meanEnergy)
#    println("Std Energy: ",stdEnergy)
    return DataSet(title,data)
end

function meanEnergy(dset::DataSet)
    return mean([i.FP_total_energy for i in dset.configs])
end

function stdEnergy(dset::DataSet)
    return std([i.FP_total_energy for i in dset.configs])
end
function rescaleData!(dataSet::DataSet,fitTo::String)
#    meanEnergy = mean([i.FP_total_energy for i in dataSet.configs])
#    stdEnergy = std([i.FP_total_energy for i in dataSet.configs])
    #offset = 3
    for i in dataSet.configs
        if fitTo == "peratom"
            i.FP_total_energy = i.FP_total_energy/ase.nAtoms(i)#(i.FP_total_energy/ase.nAtoms(i) - dataSet.meanEnergy)/dataSet.stdEnergy-dataSet.offset
       # else
       #     i.FP_total_energy = (i.FP_total_energy - dataSet.meanEnergy)/dataSet.stdEnergy-dataSet.offset
        end
    end
#    dataSet.offset = offset
#    dataSet.standardized = true
    return nothing #stdEnergy, meanEnergy
end

function writeStructuresIn(path::String, structures::DataSet)
    io = open(path, "w")
    for crystal in structures.configs
        write(io,"#-----------------------\n")
        write(io,crystal.title * "\n")
        write(io,string(crystal.latpar) * "\n")
        for lv in eachcol(crystal.lVecs)
            write(io,join(lv," ") * "\n")
        end
        write(io,join(ase.nType(crystal)," ") * "\n")
        write(io,crystal.coordSys[1] * "\n")
        for typ in crystal.atomicBasis
            for atom in typ
                write(io,join(atom," ") * "\n")
            end
        end
        write(io,"#Energy:\n")
        write(io,string(crystal.FP_total_energy) * "\n")
    end
    close(io)
    
end



end