cd("/Users/legoses/OneDrive - BYU-Idaho/codes/MatSim/src/libraries/")
push!(LOAD_PATH,pwd())

using ase
using vaspUtils
structNum = try
    parse(Int,ARGS[1])
catch y
    error("You need to specify which structure you want to build")
end
println("Building POSCAR for structure: $structNum")

currentDir = @__DIR__
path = joinpath(currentDir,"struct_enum.out")
if !isfile(path)
    error("struct_enum.out not found in current directory")
end
thisCrystal = ase.fromEnum(path,structNum,["Pt","Ag"])
vaspUtils.writePOSCAR(thisCrystal,joinpath(currentDir,"vasp." * string(structNum)))
