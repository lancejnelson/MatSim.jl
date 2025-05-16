module vaspUtils

using Printf
using ase
using DelimitedFiles




#function readStructuresIn(folder::String,file::String,species::Vector{String};overwriteLatPar = false)
#    cd(folder)
#    file = open(file,"r")
#    pos = readlines(file)
#
#    data = Vector{ase}()
#    for (idx,line) in enumerate(pos)
#
#        if occursin("#--",line)
#            nAtoms = sum([parse(Int64,x) for x in split(pos[idx + 6])])
#            startpoint = idx + 1
#            theend = idx + 7 + nAtoms
#            thisCrystal = Crystal(pos[startpoint:theend],species,overwriteLatPar = overwriteLatPar, energyPerAtomFP = parse(Float64,pos[theend + 2 ]))
#            if !isnan(thisCrystal.energyPerAtomFP)
#                push!(data,thisCrystal)
#            end
#        end
#    end
#    return DataSet(data)
#end

function readVaspFolders(folder::String,outFile::String;poscar = "CONTCAR",outcar = "OUTCAR",energy = "peratom")
    if lowercase(energy) == "fenth"
        folders = readdir(folder,join=true)
        loc = findall([occursin("pure",x) for x in folders])
        pureDirs = folders[loc]
        pures = zeros(length(pureDirs))
        if length(pureDirs) > 0
            species = sort!([String(split(x,"pure")[2]) for x in pureDirs],rev = true)
            for (idx,i) in enumerate(pureDirs)
                pure = ase.fromPOSCAR(joinpath(i,"POSCAR"),species)
                pure.title *= " (" * join(species, "-") * ")"
                pure.energyPerAtomFP = getEnergy(joinpath(i,"OUTCAR"))/pure.nAtoms
                pures[argmax(pure.nType)] = pure.energyPerAtomFP
            end
        else
            error("I can't find pure calculations so I won't be able to calculate formation energies")
        end

    end
    open(outFile, "w") do f
        write(f,energy,"\n")
        write(f,"#--------------------","\n")
    end
    for obj in readdir(folder,join = true)
        species = sort!(findSpecies(joinpath(obj,"POTCAR")),rev = true)
        if isdir(obj)
            poscarPresent = isfile(joinpath(obj,poscar))
            outcarPresent = isfile(joinpath(obj,outcar)) 
            if poscarPresent && outcarPresent
                crystal = ase.fromPOSCAR(joinpath(obj,poscar),["Pt","Ag"])
                crystal.title *= " (" * join(species, "-") * ")"
                crystal.energyPerAtomFP = getEnergy(joinpath(obj,outcar))/crystal.nAtoms
                writePOSCAR(crystal,joinpath(folder,outFile),"a")
                open(joinpath(folder,outFile), "a") do f
                    write(f,"Energy:","\n")
                    if lowercase(energy) == "peratom"
                        writedlm(f,crystal.energyPerAtomFP)
                    elseif lowercase(energy) == "total"
                         writedlm(f,crystal.energyPerAtomFP * crystal.nAtoms)
                    elseif lowercase(energy) == "fenth"
                        MatSim.formationEnergy!(crystal,pures)  # Haven't figured out how I'm gonna get pures yet.
                        writedlm(f,crystal.formationEnergyFP)
                    end
                    write(f,"#--------------------","\n")
                end
                
                
            else
                @printf("Didn't find one of the necessary files. POSCAR: %s OUTCAR: %s Skipping this dir: %s\n",poscarPresent ? "true" : "false",outcarPresent ? "true" : "false", obj)
                #print(obj)
            end
        end
    end



end


function writePOSCAR(crystal::ase.atoms,fName::String,openCode::String = "w")
    open(fName,openCode) do f
        write(f,crystal.title,'\n')
        writedlm(f,crystal.latpar',' ')
        writedlm(f,crystal.lVecs',' ')
        writedlm(f,ase.get_nTypes(crystal)',' ')
        write(f,crystal.coordSys[1],'\n')
        counter = 1
        for atom in crystal.positions
                writedlm(f,atom',' ')
               # seek(f,position(f) - 1)
            #    write(f,' ' * crystal.atomTypes[counter],'\n')
                counter +=1
        end
    
    end
    

end

function getEnergy(filePath::String)

    open(filePath) do f
        lines = readlines(f)
        energyLine = findall([occursin("free  energy",x) for x in lines])
        return parse(Float64,split(lines[energyLine[1]])[end-1])
    end

end


function findSpecies(filePath::String)

    open(filePath) do f
        lines = readlines(f)
        titleLines = findall([occursin("TITEL",x) for x in lines])
        #println(split(lines[titleLines[1]])[4])
        return [split(lines[x])[4] for x in titleLines]
    end

end

function writePOTCAR(path::String,potcar::Dict)

#    settings = YAML.load_file("VASP.yml")
    potcarsroot = potcar["path"]
    potcarspecies = reverse(sort(potcar["species"]))
    n = length(potcarspecies)
    dirs = [potcarsroot * "/" * x * "/POTCAR" for x in potcarspecies]
        
    catCommand = `cat $dirs`
    outpath = joinpath(path,"POTCAR")
    run(pipeline(catCommand, stdout = outpath))

end


function writeINCAR(writepath::String,incar::Dict)

    open(joinpath(writepath,"INCAR"), "w") do f
        for (i,d) in incar
            write(f,i * "=" * string(d) * "\n")
        end
    end
end

function writeKPOINTS(writepath::String,kpoints::Dict)

    # generate the input file.

    open(joinpath(writepath,"KPGEN"),"w") do f
        for (key,val) in kpoints["settings"][1]
            write(f,key * "=" * string(val)*"\n")
        end
        write(f,"")
    end
#    potcarsroot = potcar["path"]
#    potcarspecies = reverse(sort(potcar["species"]))
#    n = length(potcarspecies)
#    dirs = [potcarsroot * "/" * x for x in potcarspecies]
#        
    currentDir = pwd()
    cd(writepath)
    kpCommand = `kpoints.x`
    run(kpCommand)
    #wait()
    cd(currentDir)
#    outpath = joinpath(path,"POTCAR")
#    run(pipeline(catCommand, stdout = outpath))

end

#end

#myElements = ["Ni", "Pd"]
#
#atoms = [5,8]
#volume = 105
#vegardsVolume(myElements,atoms,volume)

end