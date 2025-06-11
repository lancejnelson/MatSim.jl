#!/usr/bin/env julia

cd("/Users/legoses/OneDrive - BYU-Idaho/codes/MatSim/src/libraries/")
push!(LOAD_PATH,pwd())
using Revise
#using MatSim

using StatsBase
using YAML

using enumeration
using ase
using vaspUtils
cwd = pwd()
cDir = @__DIR__

settings = YAML.load_file(joinpath(cDir,"VASP.yml"))
kp = settings["KPOINTS"]  # KPOINTs settings
dataPath = settings["path"] # Path where the data will be found
enumSettings = settings["ENUM"] # Enum settings (indicates which structures should be built)
#enum=enumeration.read_Enum_header(enumSettings["file"])  # Read the enum.out header
potcar = settings["POTCAR"]  # POTCAR settings

# Find out how many total structures are in struct_enum.out
headerLength = length(split(readuntil(joinpath(cDir,"struct_enum.out"), "start"), "\n"))
totalStructs = countlines(joinpath(cDir,"struct_enum.out")) - headerLength

# Build the list of structures to build folders for.
if enumSettings["structs"] == "random"
    structs = sample(1:totalStructs,enumSettings["nstructs"],replace=false)
elseif enumSettings["structs"] == "sequence"
    seq = parse.(Int64,split(enumSettings["nstructs"],"-"))
    structs = seq[1]:seq[2]
end

# Build each folder.
for i in structs
    println("here: ",i)
    path = joinpath(dataPath,"str" * string(i))
    mkpath(path)
#    vaspUtils.writePOTCAR(path,settings["POTCAR"])  # Write the POTCAR file
    vaspUtils.writeINCAR(path,settings["INCAR"])   # Write the INCAR file
    #enumStruct = enumeration.read_struct_from_enum(enumSettings["file"],i)
    crystal = ase.fromEnum(enumSettings["file"],i,potcar["species"])
    vaspUtils.writePOSCAR(crystal,joinpath(path,"POSCAR"),)  # Write the POSCAR file.
    vaspUtils.writeKPOINTS(path,kp)
end

