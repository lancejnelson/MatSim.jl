module ase

using StaticArrays
using LinearAlgebra
using enumeration
using LinearAlgebra:Diagonal,diag,cross,UpperTriangular,norm,det

mutable struct atoms
    title::String
    latpar::Float64
    lVecs::SMatrix{3,3,Float64,9}
#    nType:: Vector{Int64} #Number of each type of atom 
#    aType:: Vector{Int64} # Integers representing the type for each basis atom
#    nAtoms::Int64  # Total number of atoms in the unit cell
    coordSys::Vector{String} # 'D' for Direct or 'C' for cartesian
    positions:: Vector{SVector{3,Float64}}  # List of all the atomic basis vector separated by type 
    velocities:: Vector{SVector{3,Float64}} # Velocities of the atoms
    masses:: Vector{Float64}
    atomTypes:: Vector{Int64} # Integers representing the type for each basis atom
    species:: Vector{String}  # Species of the atoms 
    FP_total_energy:: Float64  # First principles total or peratom energy
    model_energy:: Float64  # Energy that should be used in the fitting process
#    formationEnergyFP:: Float64  # First-principles formation energy 
#    formationEnergyModel:: Float64  # model formation energy
#    energyPerAtomModel:: Float64 # Model total or per atom energy
 #   order::Int64 # binary, ternary, etc.
    lj_vec::Vector{Float64} # Vector of 1/r^6 and 1/r^12 values for each unique pair of atom types
#    r6::Array{Float64,2}
#    r12::Array{Float64,2}
end




latpars =Dict("H"=> 3.75,"He"=> 3.57,"Li"=> 3.49,"Be"=> 2.29,"B"=> 8.73,"C"=> 3.57,"N"=> 4.039,
              "O"=> 6.83,"Ne"=> 4.43,"Na"=> 4.23,"Mg"=> 3.21,"Al"=> 4.05,"Si"=> 5.43,"P"=> 7.17,
              "S"=> 10.47,"Cl"=> 6.24,"Ar"=> 5.26,"K"=> 5.23,"Ca"=> 5.58,"Sc"=> 3.31,"Ti"=> 2.95,
              "V"=> 3.02,"Cr"=> 2.88,"Mn"=> 8.89,"Fe"=> 2.87,"Co"=> 2.51,"Ni"=> 3.52,"Cu"=> 3.61,
              "Zn"=> 2.66,"Ga"=> 4.51,"Ge"=> 5.66,"As"=> 4.13,"Se"=> 4.36,"Br"=> 6.67,"Kr"=> 5.72,
              "Rb"=> 5.59,"Sr"=> 6.08,"Y"=> 3.65,"Zr"=> 3.23,"Nb"=> 3.3,"Mo"=> 3.15,"Tc"=> 2.74,
              "Ru"=> 2.7,"Rh"=> 3.8,"Pd"=> 3.89,"Ag"=> 4.09,"Cd"=> 2.98,"In"=> 4.59,"Sn"=> 5.82,
              "Sb"=> 4.51,"Te"=> 4.45,"I"=> 7.27,"Xe"=> 6.2,"Cs"=> 6.05,"Ba"=> 5.02,"Hf"=> 3.2,
              "Ta"=> 3.31,"W"=> 3.16,"Re"=> 2.76,"Os"=> 2.64,"Ir"=> 3.84,"Pt"=> 3.92,"Au"=> 4.08,
              "Hg"=> 2.99,"Tl"=> 3.46,"Pb"=> 4.95,"Bi"=> 4.75)

element_volume =Dict("H"=>37.2958,"He"=>32.1789,"Li"=>21.2543,"Be"=>8.49323,"B"=>7.24205,"C"=>5.68741,
                 "N"=>46.6002,"O"=>22.2802,"F"=>17.0258,"Ne"=>21.7346,"Na"=>23.2596,"Mg"=>23.3928,
                 "Al"=>16.6075,"Si"=>7.8511,"P"=>9.1459,"S"=>17.1672,"Cl"=>35.2074,"Ar"=>36.3829,
                 "K"=>71.5278,"Ca"=>43.4353,"Sc"=>25.6478,"Ti"=>18.1565,"V"=>13.7718,"Cr"=>11.9439,
                 "Mn"=>19.3207,"Fe"=>11.82,"Co"=>11.1838,"Ni"=>10.9036,"Cu"=>11.7615,"Zn"=>13.311,
                 "Ga"=>18.4496,"Ge"=>19.3638,"As"=>20.4270,"Se"=>58.6173,"Br"=>33.3170,"Kr"=>46.7873,
                 "Rb"=>87.3384,"Sr"=>56.1889,"Y"=>33.0792,"Zr"=>23.8327,"Nb"=>17.9685,"Mo"=>15.6279,
                 "Tc"=>14.5458,"Ru"=>13.9206,"Rh"=>13.718,"Pd"=>14.716,"Ag"=>17.1045,"Cd"=>18.7161,
                 "In"=>26.6861,"Sn"=>29.3238,"Sb"=>27.1733,"Te"=>62.3227,"I"=>24.3807,"Xe"=>59.582,
                 "Cs"=>110.723,"Ba"=>63.253,"Hf"=>23.1748,"Ta"=>18.1323,"W"=>15.7772,"Re"=>14.8694,
                 "Os"=>14.5485,"Ir"=>14.1558,"Pt"=>15.0591,"Au"=>16.9793,"Hg"=>27.6914,"Tl"=>29.2949,
                 "Pd"=>30.3218,"Bi"=>31.2849,"U"=>13.6389, "O"=>13.6389)#"U"=>50.84, "O"=>15.86}


#config(filePath::String, species::Vector{String}; overwriteLatPar = false) = fromPOSCAR(filePath,species,overwriteLatPar = overwriteLatPar)
#config(list::Vector{String}, species::Vector{String}; overwriteLatPar = false) = fromPOSCAR(list, species,overwriteLatPar = overwriteLatPar)
#config(filePath::String, strN::Int64,species::Vector{String};mink::Bool=true) = fromEnum(enum,enumStruct,species,mink=mink)


function index_to_integer(iType,jType,order)
    a = minimum([iType,jType])
    b = maximum([iType,jType])
    return (div((a-1) * (2 * order + 2 - a),2)) + (b - a) + 1
end

function integer_to_index(k::Int, n::Int)::Union{Tuple{Int,Int},Nothing}
    # Check if k is valid
    max_k = (n * (n + 1)) ÷ 2
    if k < 1 || k > max_k
        return nothing
    end

    # Find a such that k falls in the range for min(i,j) = a
    for a in 1:n
        # Pairs up to i' = a-1
        s = ((a - 1) * (2 * n + 2 - a)) ÷ 2
        # Number of pairs for i' = a (j' = a to n)
        count = n - a + 1
        # Check if k is in range [s + 1, s + count]
        if s < k ≤ s + count
            # Compute b
            b = k - s + a - 1
            return (a, b)
        end
    end

    return nothing  # Should not reach here if k is valid
end


function nAtoms(atoms::atoms)
    """
    Returns the total number of atoms in the crystal.
    """
    return sum([length(a) for a in atoms.positions])
end

function nType(atoms::atoms)
    """
    Returns the number of different types of atoms in the crystal.
    """
    return [length(a) for a in atoms.positions]
end

function rescaleEnergy(energy,mean_energy, std_energy,offset)

    return (energy - meanEnergy)/std_energy - offset

end

function undo_rescaleEnergy(energy,mean_energy, std_energy,offset)
    return (energy + offset) * std_energy + mean_energy

end
"""
    fromEnum(enum, enumStruct,species, [, energyPerAtomFP, modelEnergy, mink])

    Builds a Crystal from the general enumeration data (enum) and the specific enumeration data for a given superstructure (enumStruct)
    `enum` contains information about the parent lattice, whereas `enumStruct` contains all of the information about a given superstructure (HNF, labeling, etc.)
    
    In addition, you can specify an energy to attach to the crystal if you choose and you can specify if you want the parent lattice vectors to be "minkowski'd"
"""
#function fromEnum(enum::Enum,enumStruct::EnumStruct,species:: Vector{String};mink=true)
function fromEnum(file::String,strN::Int64,species:: Vector{String};mink=true)
    enum = enumeration.read_header(file) 
    enumStruct = enumeration.read_struct(file,strN)

    cardinalDirections = Float64[0 0 0
                          1 0 0
                          -1 0 0
                          0 1 0
                          0 -1 0
                          0 0 1
                          0 0 -1]
    sLV = enum.pLV * enumStruct.HNF  # Get super lattice vectors!
    a,b,c,d,e,f = enumStruct.HNF[1,1],enumStruct.HNF[2,1],enumStruct.HNF[2,2],enumStruct.HNF[3,1],enumStruct.HNF[3,2],enumStruct.HNF[3,3]
    # Get atomic basis..
    nAtoms = length(enumStruct.labeling)
    aBas = [zeros(3) for i =1:nAtoms]
    gIndex = zeros(Int64,nAtoms)
    iC = 0
    for iD = 1:enum.nD
        for z1 = 0:a-1, z2= Int((b * z1) / a):Int(c + (b * z1) / a) - 1, z3 = Int(z1 * (d - (e * b) / c) / a + (e * z2) / c):Int(f + z1 * (d - (e * b) / c) / a + (e * z2) / c) - 1
            iC += 1
            temp = enum.pLV * [z1,z2,z3]
            temp2 = temp + enum.dVecs[iD]
            aBas[iC] .= temp2

            gReal = enumStruct.L * [z1,z2,z3]

            g = Int.(enumStruct.L * [z1,z2,z3])
            if !(g ≈ gReal)
                error("Mapping didn't work")
            end

            # Bring g-vector back into the first cell
            g = g .% diag(enumStruct.SNF)
            # Which character in the labeling string corresponds to this atom..
            gIndex[iC] = Int((iD - 1) * enumStruct.SNF[1,1] * enumStruct.SNF[2,2] * enumStruct.SNF[3,3]
                         + g[1] * enumStruct.SNF[2,2] * enumStruct.SNF[3,3]
                         + g[2] * enumStruct.SNF[3,3]
                         + g[3]) + 1
            
        end
    end
    aType = [parse(Int,enumStruct.labeling[i]) for i in gIndex]  # Get the atomic labels in the right order (i.e. using gIndex)
    nType = [count(==(i),aType) for i in sort!(collect(Set(aType)))]
    if !all(isapprox.(sort(gIndex), gIndex,atol = 1e-7))
        println(aType)
        println(sort(aType))
        println(gIndex)
        #println("Found one")
        error("This isn't an error. I was just curious if the gIndex Vector is ever in an odd order.")
    end
    idx = sortperm(aType)  # Sort the type list, saving the indices
    #aBas = getPartitions(aBas[idx],nType) #Sort the atomic basis to match and then partition them by atom type.
    aType = sort(aType) .+ 1

    # The next two lines are to add the displacements in using the arrows string from struct_enum.out. It's not tested
    # yet and I need to verify with Gus that I'm doing it right.  Just delete the next three lines to undo it if you find that
    # it's not right.
    arrows = [parse(Integer,enumStruct.arrows[i]+1) for i in gIndex]  #Get the list of arrows (as integers) in the right order.
    displacements = [cardinalDirections[i,:] for i in arrows] # Partition to the same shape as the basis list so we can easily add them.
    aBas .+= 0.1 * displacements

    cellVolume = abs(cross(sLV[:,1],sLV[:,2])' * sLV[:,3])
    
    nInteractions = Int(enum.k * (enum.k + 1)/2)

#    r6 = UpperTriangular(zeros(enum.k,enum.k))
#    r12 = UpperTriangular(zeros(enum.k,enum.k))
    lj_vec = zeros(2*nInteractions)
    masses = [0.0 for i in 1:nAtoms]
    velocities = [SVector{3,Float64}(zeros(3)) for i in 1:nAtoms]
    final = atoms(enum.title * " str #: " * string(enumStruct.strN),1.0,mink ? minkowski_reduce(sLV,1e-5) : sLV,["C"],aBas,velocities, masses,aType,["Unk" for i in nType],0.0,0.0,lj_vec)
    CartesianToDirect!(final)
    final.latpar = vegardsVolume(species,nType,cellVolume)
    return final
end

"""
When building crystals we often start by building a list of atomic coordinates, irrespective of which type of atom they are.
When you're all done, you'd like to partition that list by atomic type.  This is handy for accessing them later because the first index
is the atom type and the second index is the number in that sublist.  This function simply performs that partitioning.

"""
function getPartitions(atoms,nAtoms::Vector{Int64})#::Vector{Vector{Float64}}
    if length(atoms) != sum(nAtoms)
        println("Number of atoms doesn't match")
        return
    end
    parts = prepend!(cumsum(nAtoms),0)
    return [atoms[parts[i]+1:parts[i+1]] for i=1:length(parts)-1]
end

"""
The starting point for a nested sampling algorithm is a set of random configurations.  This function builds a single atomic configuration
with atomic positions at random locations.

# Arguments
- `lPar: Float64`: lattice parameter
- `lVecs: Matrix{Float64}`: lattice vectors (by columns)
- `nAtoms: Vector{Int64}`: Vector containing the number of each atom type in the configuration
- `cutoff: Float64`: No two atoms can be closer than this value.
- `species: Vector{String}`: Vector of Strings representing the atomic species.  Might not be necessary!


"""
# This function initializes only the cell parameters for an atoms object. All other type parameters are set to zero.
function initialize_cell_shape(lPar::Float64, lVecs::SMatrix{3,3,Float64,9},nAtoms::Int)

    return atoms("No Atoms", lPar,lVecs,["C"],[SVector{3,Float64}(zeros(3)) for x= 1:nAtoms],[SVector{3,Float64}(zeros(3)) for x = 1:nAtoms],[1.0 for x=1:nAtoms],[1 for x = 1:nAtoms],["Unk"],0.0,0.0,zeros(3))

end

function set_coord_sys!(atoms,coordsys)
    #println(lowercase(coordsys)[1])
    if lowercase(coordsys)[1] == 'c'
        #println("Here1")
        DirectToCartesian!(atoms)
    elseif lowercase(coordsys)[1] == 'd'
        #println("here2")
        CartesianToDirect!(atoms)
    end
end 

function set_cell!(atoms,new_cell; scale_atoms = false)
    # Note: This routine is not meant to preserve the periodicity of the original configuration.  In other words, the number of atoms inside the cell will not change no matter how much the volume changes.  

    oldCell = atoms.lVecs
    atoms.lVecs = new_cell  # Can we find a way to assign in-place.  Can't be an SMatrix if so.
    currentCoordSys = atoms.coordSys[1]
    T = new_cell * inv(atoms.lVecs)  # Transformation matrix

    if scale_atoms
        DirectToCartesian!(atoms)
        for (iAtom,atom) in enumerate(atoms.positions)
            atoms.positions[iAtom] =  T * atoms.positions[iAtom]
        end
    end
    set_coord_sys!(atoms,currentCoordSys)

end

function cell_volume(cell)
    return abs(cell[:,1]' * cross(cell[:,2],cell[:,3]))
end

function min_aspect_ratio(cell)
    volume = cell_volume(cell)
    min_aspect = 1000
    for i=1:3
        vi = cell[i,:]
        vnorm_hat = cross(cell[i%3 + 1,:],cell[(i+1)%3 + 1,:])
        vnorm_hat /= norm(vnorm_hat)
        min_aspect = minimum([min_aspect,abs(vnorm_hat' * vi)])

    end
    min_aspect /= volume^(1.0/3.0)
    return min_aspect

end


function get_concentrations(atoms::atoms)
    """
    Returns the concentrations of each atom type in the crystal.
    """
    order = length(unique(atoms.atomTypes))
    numAtoms = nAtoms(atoms)
    concentrations = [0.0 for i=1:order]
    for i=1:order
        concentrations[i] = length(atoms.positions[i])/numAtoms
    end
    return concentrations
end


# This function will take an already created atoms object and randomly place atoms in the cell.
#function set_atoms_random!(atoms::atoms,nAtoms::Vector{Int64},cutoff::Float64,species::Vector{String})
#    initialize_with_random_positions!(atoms,nAtoms,cutoff,species)
    #atoms.nType = newAtoms.nType
    #atoms.atomTypes = newAtoms.aType
    #atoms.nAtoms = newAtoms.nAtoms
    #atoms.coordSys = newAtoms.coordSys
    #atoms.positions = newAtoms.positions
    #atoms.species = newAtoms.species
    #atoms.order = newAtoms.order
    #atoms.
    #atoms.title = join(species,"-")
#end

# This routine will initialize an atoms object so that the atoms are randomly placed in the cell.
function set_atoms_random!(atoms,nAtoms::Vector{Int64},min_separation::Float64, cutoff,species::Vector{String})

    lPar = atoms.latpar
    lVecs = atoms.lVecs
    totalAtoms = sum(nAtoms)
    order = length(nAtoms)
    aTypes = hcat([ [n for i=1:nAtoms[n]]' for n=1:length(nAtoms)]...)'
    loopBounds = SVector{3,Int64}(convert.(Int64,cld.(cutoff ,SVector{3,Float64}(norm(x) for x in eachcol(atoms.latpar * atoms.lVecs)) )))
    atoms.positions = []
    #atoms.positions .= [zeros(Float64,3) for i=1:totalAtoms]  
    nTotal = 0
    counter = 0
    while length(atoms.positions) < totalAtoms
        counter += 1
        if counter > 5000
            println("Having trouble putting that many atoms into this simulation cell.")
            println("So far I've only place $(length(atoms.positions)) atoms." )
            println(counter)
            error("Stopping")
        end
        newAtomCart = lPar * lVecs * rand(3)
        println("new Cart", newAtomCart)
        if !atom_inside_cutoff(atoms,newAtomCart,min_separation)
            push!(atoms.positions, CartesianToDirect(atoms.latpar * atoms.lVecs,SVector{3,Float64}(newAtomCart)))
        end
        display(atoms.positions)
        continue
        #println("considering new atom: " ,newAtomCart)
        #println("All previous atoms:")
        #display(positions[1:nTotal])
        newOK = true
        for i=1:nTotal
            println("checkin min distances")
            display(atoms.positions)
            println(i)
            println(get_single_atom_neighbors(atoms,cutoff,atoms.positions[i],loopBounds))
            println("n Total", nTotal)
            #if !atom_inside_cutoff(atoms,newAtomCart,cutoff)
            if min_distance < min_separation
                newOK = false
            end
            if !newOK
          #      println("rejecting at index: ", i)
                break
            end
        end
        if newOK
         #   println("accepting")
            nTotal += 1
            push!(atoms.positions, SVector{3,Float64}(newAtomCart))
#            atoms.positions[nTotal] = newAtomCart
        end
    end
    #positions = getPartitions(atom_locations,nAtoms)
    k = length(species)
#    r6 = UpperTriangular(zeros(k,k))
#    r12 = UpperTriangular(zeros(k,k))
    lj_vec = zeros(2*Int(k*(k+1)/2))
    #println("lj_vec: ",lj_vec)
    atoms.masses = [0.0 for i in 1:length(atoms.positions)]
    atoms.velocities = [SVector{3,Float64}(zeros(3)) for i in 1:length(atoms.positions)]
    atoms.species = species
    atoms.coordSys = ["D"]
    atoms.lj_vec = lj_vec
    atoms.title = "Random Locations"
    atoms.atomTypes = aTypes
    #    return atoms("Random Locations",lPar,lVecs,["C"], positions,velocities,masses,species,0.0,0.0,lj_vec)
end

"""
Construct a Crystal from a POSCAR file.

# Arguments
- `filePath:: String`: Location of POSCAR
- `species:: Vector{String}`: Vector of Strings representing the atomic species.  Might not be necessary!
- `overwriteLatPar:: Bool = false`: Do you want to keep the latpar in the file or compute it using Vegard's law

"""

function fromPOSCAR(filePath::String,species::Vector{String};overwriteLatPar = false)

#    cd(folder)
    
    file = open(filePath, "r")
    pos = readlines(file)

    title = pos[1]
    lVecs = SMatrix{3,3}(reduce(hcat,[parse.(Float64, y) for y in [split(x) for x in pos[3:5]]])) # Read in the lattice vectors.
    if !isletter(lstrip(pos[6])[1])
        counter = 6
    else
        counter = 7
    end
    nBasis = parse.(Int,split(pos[counter]))
    coordSys = [pos[counter + 1]]
    latpar = parse(Float64,pos[2])
    aType = hcat([ [n for i=1:nBasis[n]]' for n=1:length(nBasis)]...)'
    positions = [parse.(Float64,split(x)[1:3]) for x in pos[(counter + 2):(counter + 1 +sum(nBasis))]] # Read all of the basis vectors SVector{3,Float64}(
    allTypes = try
        [split(x)[end] for x in pos[(counter + 2):end]]
    catch e
        println("No atom types listed in the POSCAR")
        ["?" for x in pos[(counter + 2):end]]
    end
#    positions = getPartitions(allBasis,nBasis)
    order = length(nBasis)
    nAtoms = sum(nBasis)


    nInteractions = Int(order * (order + 1)/2)
    # I'm not sure why I'm putting this code here. We'll always want to use
    # the lattice parameter found in the structures.in file because that
    # was the geometry used in the VASP calculation. 
    cellVolume = abs(cross(lVecs[:,1],lVecs[:,2])' * lVecs[:,3])
    if overwriteLatPar
        println("Overwriting lat par in file with Vegard's law")
        latpar = vegardsVolume(species,nBasis,cellVolume)
    else
#        println("Keeping lattice parameter in file")
    end
    
#    r6 = UpperTriangular(zeros(order,order))
#    r12 = UpperTriangular(zeros(order,order))
    lj_vec = zeros(2*nInteractions)
#    fEnth = formationEnergy(pureEnergies,nBasis ./ nAtoms,energyPerAtomFP)
    masses = [0.0 for i in 1:nAtoms]
    velocities = [SVector{3,Float64}(zeros(3)) for i in 1:nAtoms]

    return atoms(title, latpar,lVecs,coordSys,positions,velocities,masses,aType,species,0.0,0.0,lj_vec)  # Create new crystal object.

end


"""
Construct a Crystal from a list of lines (like when reading a structures.in file and you pick off the POSCARs one by one.)

# Arguments
- `lines:: Vector{String}`: POSCAR lines
- `species:: Vector{String}`: Vector of Strings representing the atomic species.  Might not be necessary!
- `overwriteLatPar:: Bool = false`: Do you want to keep the latpar in the file or compute it using Vegard's law

"""

function fromPOSCAR(lines::Vector{String},species::Vector{String};overwriteLatPar = false)
    title = lines[1]
    lVecs = reduce(hcat,[parse.(Float64, y) for y in [split(x) for x in lines[3:5]]])#SMatrix{3,3}() # Read in the lattice vectors.
    if !isletter(lstrip(lines[6])[1])
        counter = 6
    else
        counter = 7
    end
    nBasis = parse.(Int,split(lines[counter]))
    coordSys = [lines[counter + 1]]
    latpar = parse(Float64,lines[2])
    aType = hcat([ [n for i=1:nBasis[n]]' for n=1:length(nBasis)]...)'
    positions = [parse.(Float64,split(x)[1:3]) for x in lines[(counter + 2):(counter + 1)+sum(nBasis)]] # Read all of the basis vectors
    allTypes = try
        [split(x)[end] for x in lines[(counter + 2):end]]
    catch e
        println("No atom types listed in the POSCAR")
        ["?" for x in pos[(counter + 2):end]]
    end
#    positions = getPartitions(allBasis,nBasis)
    order = length(nBasis)
    nAtoms = sum(nBasis)

    nInteractions = Int(order * (order + 1)/2)
    # I'm not sure why I'm putting this code here. We'll always want to use
    # the lattice parameter found in the structures.in file because that
    # was the geometry used in the VASP calculation. 
    cellVolume = abs(cross(lVecs[:,1],lVecs[:,2])' * lVecs[:,3])
    if overwriteLatPar
        println("Overwriting lat par in file with Vegard's law")
        latpar = vegardsVolume(species,nBasis,cellVolume)
    else
#        println("Keeping lattice parameter in file")
    end
#    r6 = UpperTriangular(zeros(order,order))
#    r12 = UpperTriangular(zeros(order,order))
    lj_vec = zeros(2*nInteractions)

    masses = [0.0 for i in 1:nAtoms]
    velocities = [SVector{3,Float64}(zeros(3)) for i in 1:nAtoms]
#    fEnth = formationEnergy(pureEnergies,nBasis ./ nAtoms,energyPerAtomFP)
    return atoms(title, latpar,lVecs,coordSys,positions,velocities,masses,aType,species,0.0,0.0,lj_vec)  # Create new crystal object.

end

function mapIntoCell!(atoms::atoms)
    CartesianToDirect!(atoms)

    # Is there a better way to do this than using the new_point variable and then
    # turning it into an SVector at the end. I can't modify the SVectors in place
    new_point = [0.0,0.0,0.0]
    for (iAtom,atom) in enumerate(atoms.positions)
        #for (iAtom,atom) in enumerate(atomType)
        for (iCoord,coord) in enumerate(atom)
            if coord < 0.0 || coord > 1.0
            #   print(i,' i')
            #   print(floor(i),' floor')
            #   print(i - floor(i),' result')
                new_point[iCoord] = coord - floor(coord)
                #atoms.positions[iType][iAtom][iCoord] = coord - floor(coord)
            elseif isapprox(coord,1.0,atol = 1e-4)#  i == 1.0
                new_point[iCoord] = 0.0
               # atoms.positions[iType][iAtom][iCoord] = 0.0
            else
                new_point[iCoord] = coord
            end
        end
        atoms.positions[iAtom] = SVector{3,Float64}(new_point)
    end

    DirectToCartesian!(atoms)
end


"""
     Calculate the formation energy for a crystal

# Arguments
- `pureEnergies:: Vector{Float64}`: Energies of pure crystals
- `concs:: Vector{Float64}`: Vector of concentrations for the crystal being considered
- `energyPerAtom:: Float64`: Energy/atom of the crystal

"""

function fccPures(types)

    lVecs = @SMatrix [0.5 0.5 0
             0.5 0 0.5
             0 0.5 0.5]
    lP = [latpars[x] for x in sort!(types,rev = true)]
    nType = [1, 0]
    aType = [1, 0]
    nAtoms = 1
    coordSys = ["D"]
    positions = [@SVector [0., 0 , 0]]
    species = sort!(types,rev =true)
    order = length(types)
#    r6 = UpperTriangular(zeros(order,order))
#    r12 = UpperTriangular(zeros(order,order))
    nInteractions = Int(order * (order + 1)/2)
    lj_vec = zeros(2 * nInteractions)
    title = join(types, "-")
    println(title)
    masses = [0.0 for i in 1:nAtoms]
    velocities = [SVector{3,Float64}(zeros(3)) for i in 1:nAtoms]
    return atoms(title,lP[1],lVecs,coordSys,positions,velocities,masses,aType,species,0.0,0.0,lj_vec),
           atoms(title,lP[2],lVecs,coordSys,positions,velocities,masses,aType,species,0.0,0.0,lj_vec) 

end

function formationEnergy(mixEnergyPerAtom,pureEnergiesPerAtom,concentrations)
   #println("Calculating formation energy")
   #println(mixEnergyPerAtom)
   #println(pureEnergiesPerAtom)
   #println(concentrations)
    return mixEnergyPerAtom - sum(concentrations .* pureEnergiesPerAtom)

end

#function formationEnergyModel!(crystal,pures)
#     crystal.formationEnergyModel = crystal.modelEnergy - sum(crystal.nType /crystal.nAtoms .* pures)
#     crystal.formationEnergyFP = crystal.energyPerAtomFP - sum(crystal.nType /crystal.nAtoms .* pures)
##    return energyPerAtom - sum(concs .* pureEnergies)
#
#end
#
#function formationEnergyFP!(crystal,pures)
#    crystal.formationEnergyFP = crystal.energyPerAtomFP - sum(crystal.nType /crystal.nAtoms .* pures)
##    return energyPerAtom - sum(concs .* pureEnergies)
#
#end

function totalEnergyFromFormationEnergy!(atoms,pures)
#    crystal.formationEnergyFP = crystal.energyPerAtomFP - sum(crystal.nTypes/crystal.nAtoms .* pures)
    atoms.energyPerAtomFP = atoms.formationEnergyFP + sum(atoms.nTypes/atoms.nAtoms .* pures)
#    return energyPerAtom - sum(concs .* pureEnergies)

end

function DirectToCartesian!(atoms::atoms)
#    println(atoms.coordSys[1])
    if lowercase(atoms.coordSys[1])[1] == 'd'
        #println("Converting to cartesian")
        atoms.positions .= [atoms.latpar * atoms.lVecs * i for i in atoms.positions ]
        #println(typeof(crystal.positions[1][1]))

        atoms.coordSys[1] = "Cart"
 #   else
 #       println("Already in Cartesian coordinates")
    end
end

function cell_volume(atoms::atoms)
    return atoms.latpar^3 * abs(cross(atoms.lVecs[:,1],atoms.lVecs[:,2])' * atoms.lVecs[:,3])
end

function CartesianToDirect!(atoms::atoms)
    if lowercase(atoms.coordSys[1])[1] == 'c'
        #println("Converting to direct")
        atoms.positions .= [ round.( ( inv(atoms.latpar * atoms.lVecs) * i) .% 1,sigdigits = 8) for i in atoms.positions ]
        #println(typeof(atoms.positions[1][1]))
        atoms.coordSys[1] = "Direct"
    #else
       # println("Already in Cartesian coordinates")
    end
end

function DirectToCartesian(lVecs::SMatrix{3,3,Float64},atom:: SVector{3,Float64})
    return lVecs * atom
end

function CartesianToDirect(lVecs::SMatrix{3,3,Float64},atom:: SVector{3,Float64})
    return inv(lVecs) * atom  .% 1
end

function isEquivalentAtom(atomOne,atomTwo)  # Atoms are assumed to be in direct coordinates.

    #=If the atoms are separated by a lattice vector we can be sure that they are equivalent atoms (translationally).
       To check this, subtract the two vectors and see if the resulting vector is only composed of 
       integers. 
    =#
    return all(round.(Int64,atomOne - atomTwo) - (atomOne - atomTwo) == [0,0,0]) 


end




function Gaussian_Reduce(U,V,eps)  
    """This routine takes two vectors (in three-space) and reduces them to
    form a shortest set (Minkowski reduced). The idea is to subtract
    multiples of U from V so that the new V is as close to the origin
    as any lattice point along the line that passes through U in the
    direction of V. The process is repeated until the new vector isn't
    shorter than the other. It's pretty obvious if you do an example
    by hand. Also see 3.1 of Lecture notes in computer science, ISSN
    0302-974, ANTS - VI: algorithmic number theory, 2004, vol. 3076,
    pp. 338-357 ISBN 3-540-22156-5. Fixes made Apr 2012 GLWH (not sure
    if they made a practical difference though)
    :arg U: a vector
    :arg V: a vector
    :arg eps: finite precision tolerance
    """

#    from numpy.linalg import norm
#    from numpy import dot

    it = 0
    if norm(U) > (norm(V) - eps)
        # Make sure that the {U,V} are listed in ascending order; ||U||<||V||
        temp = copy(U)
        U = copy(V)
        V = copy(temp)  # Keep V as the longest vector
    end
    done = false
    it = 1
    while  !done
        if it > 10  # pragma: no cover
            error("gaussian_reduce_two_vectors failed to converge in 10 iterations")
        end
        R = V - Int(round((U'*V)/(U'*U) + 1e-10)) * U

        V = copy(U)  # Swap U and V (so U remains the shortest)
        U = copy(R)
        if norm(U) >= (norm(V) - eps)
            done = true
        end
        it += 1
    end
    # Make sure that the {U,V} are listed in ascending order on exit; ||U||<||V||
    temp = copy(U)
    U = copy(V)
    V = copy(temp)
    return U, V
end


function reduce_C_in_ABC(ABC,eps)
    #oldABC = deepcopy(ABC)
    A,B,C = ABC[:,1],ABC[:,2], ABC[:,3]
    A,B = Gaussian_Reduce(A,B,eps)
    ABC = [A B C]  # Update ABC to have the updated vectors in it.
    cpdAB = cross(A,B)/norm(cross(A,B))  # Find the unit vector the points perp to A-B plane

    T = C - cpdAB * (C' * cpdAB)  # Find projection of C onto A-B plane
    if !isapprox(T' * cross(A,B),0,atol = eps)  # Check to see if we did it right. Dot product should be zero now.
        error("Projection of C into A,B plane failed.")
    end

#    ABC = hcat(A, B, C)  # lVs by columns, not rows.
    LC = Int.(floor.(inv(ABC) * T .+ eps)) # Find what combinations of A,B,C will produce the projection.

    corners = [0 0 0
               1 0 0
               0 1 0
               1 1 0]  # The projection is closest to one of these corners.
    distances = [norm(T - ABC * (corners[i,:] + LC) ) for i =1:4]  # Get the distances to each of the corners, so I know which combination of ABC to subtract off.
    val,idx = findmin(distances)

    C = C - (ABC * (corners[idx,:] + LC))  # This is guaranteed to produce another lattice vector because you are subtracting multiples of other lattice vectors.

    newABC = hcat(A,B,C)
    if !all(isapprox.(inv(newABC) * ABC - Int.(round.(inv(newABC) * ABC)),0,atol = eps))
        display(inv(newABC)* ABC)
        display(Int.(round.(inv(newABC) * ABC)))
        println("This matrix should have been integers")
        error("Lattice not preserved in reduce_C_in_ABC")
    end
    return newABC

end





function get_single_atom_force_lj(model,atoms::atoms,centerAtom::Int64, loopBounds::SVector{3,Int64})
    #ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to k-nary case.
    CartesianToDirect!(atoms)
    order = length(unique(atoms.atomTypes))

    addVec = zeros(3)
    indices = zeros(2)
    fVec = SVector(0,0,0)
    centerType = atoms.atomTypes[centerAtom]
    for (iAtom,neighboratom) in enumerate(atoms.positions)  #Loop over the different atom types.
        neighborType = atoms.atomTypes[iAtom]
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            addVec .= (i,j,k) 
            newAtom = neighboratom + addVec  #Periodic image of this atom
            newCart = DirectToCartesian(atoms.latpar * atoms.lVecs,newAtom)  # Convert to cartesian coordinate system
            r = newCart - DirectToCartesian(atoms.latpar * atoms.lVecs,atoms.positions[centerAtom]) 
            if norm(r) < model.cutoff && !isapprox(norm(r),0,atol = 1e-3)
                #println("Adding to force")
#                indices = iNeighbor < centerAtom[1] ? @SVector[iNeighbor,centerAtom] : @SVector[centerAtom,iNeighbor]
                index = index_to_integer(centerType,neighborType,order)
                fVec -=    12. * model.ϵ[index] * model.σ[index]^12/norm(r)^13 * r/norm(r)
                fVec +=    6. *  model.ϵ[index] * model.σ[index]^6/norm(r)^7 * r/norm(r)
            end
        end
    end
    return fVec
end


function get_single_atom_LJ_PE(model,atoms::atoms,centerAtom::SVector{3,Float64}, centerType:: Integer, loopBounds::SVector{3,Int64})
    get_single_atom_r6_r12!(atoms,model.cutoff,centerAtom,centerType,loopBounds)
    coeffs = vcat(-1.0 * model.ϵ .* model.σ.^6, model.ϵ .* model.σ.^12)
    totalEnergy = atoms.lj_vec' * coeffs
    
    return (totalEnergy + model.offset) * model.stdEnergy + model.meanEnergy
end

function get_neighbor_distances(atoms::atoms, cutoff)
    CartesianToDirect!(atoms)
    
    eVals = eigvals(transpose(atoms.latpar .* atoms.lVecs) * (atoms.latpar .* atoms.lVecs))
    maxN = cutoff/sqrt(minimum(eVals))
    
    loopBounds = SVector{3,Int64}(convert.(Int64,cld.(cutoff ,SVector{3,Float64}(norm(x) for x in eachcol(atoms.latpar * atoms.lVecs)) )))
    # The outer two loops are to loop over different centering atoms.
    nn_distances = []
    for (iCenter,centerAtom) in enumerate(atoms.positions)
        centerAtomC = DirectToCartesian(atoms.latpar * atoms.lVecs,centerAtom)
        push!(nn_distances,get_single_atom_neighbors(atoms,cutoff,centerAtomC,loopBounds))    # Find the contribution to the LJ energy for this centering atom.
        
    end

    return nn_distances
end

function atom_inside_cutoff(atoms::atoms,centerAtom::SVector{3,Float64},cutoff)
    loopBounds = SVector{3,Int64}(convert.(Int64,cld.(cutoff ,SVector{3,Float64}(norm(x) for x in eachcol(atoms.latpar * atoms.lVecs)) )))
    println(loopBounds)
    for neighboratom in atoms.positions  #Loop over the different atom types.
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            newAtom = neighboratom + @SVector[i,j,k]
            newCart = DirectToCartesian(atoms.latpar * atoms.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < cutoff && !isapprox(r,0.0,atol = 1e-3)
                return true
            end
        end
    end
    return false
end

function get_single_atom_neighbors(atoms::atoms,cutoff,centerAtom::SVector{3,Float64}, loopBounds::SVector{3,Int64})
    nns = []
    for (iNeighbor,neighboratom) in enumerate(atoms.positions)  #Loop over the different atom types.
        neighborType = atoms.atomTypes[iNeighbor]
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            newAtom = neighboratom + @SVector[i,j,k]
            newCart = DirectToCartesian(atoms.latpar * atoms.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < cutoff && !isapprox(r,0.0,atol = 1e-3)
                push!(nns,(iNeighbor,r))
            end

        end
    end

    return nns
end


function get_single_atom_r6_r12!(atoms::atoms,cutoff,centerAtom::SVector{3,Float64}, centerType:: Integer, loopBounds::SVector{3,Int64})
    #    ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    order = length(unique(atoms.atomTypes))
    nInteractionTypes = sum(1:order)
    foundInteractions = 0
    for (iNeighbor,neighboratom) in enumerate(atoms.positions)  #Loop over the different atom types.
        neighborType = atoms.atomTypes[iNeighbor]
            # And these three inner loops are to find all of the periodic images of a neighboring atom.
        for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]
            newAtom = neighboratom + @SVector[i,j,k]
            newCart = DirectToCartesian(atoms.latpar * atoms.lVecs,newAtom)
            r = norm(newCart - centerAtom) 
            if r < cutoff && !isapprox(r,0.0,atol = 1e-3)
                foundInteractions += 1
                #If the neighbor atom is inside the unit cell, then its going to be
                # double counted at some point when we center on the other atom.  
                # So we count it as half each time.
                # The LJ parameters are stored in the upper triangular portion of a matrix
                # The bottom triangle is redundant.. interaction between a and b is equivalent
                # to interaction between b and a.  So I sort the indices here so that the bottom
                # triangle of the matrix never gets updated, only the upper right.
                #indices[1],indices[2] = iNeighbor,centerType
                #println(atoms.atomTypes)
                index = index_to_integer(centerType,neighborType,order)
                #println(centerType,neighborType,order,index)
                if all(@SVector[i,j,k] .== 0 )
                    atoms.lj_vec[index] +=  4.0 * 1.0/2.0 * 1.0/r^6
                    atoms.lj_vec[index + nInteractionTypes] +=   4.0 * 1.0/2.0 * 1.0/r^12
                else
                    atoms.lj_vec[index] +=  4.0 * 1.0/r^6
                    atoms.lj_vec[index + nInteractionTypes] +=  4.0 *  1.0/r^12
                end
            #    indices = iNeighbor < centerType ? @SVector[iNeighbor,centerType] : @SVector[centerType,iNeighbor]
            #    if all(@SVector[i,j,k] .== 0 ) 
            #        atoms.r6[indices[1],indices[2]] +=  4. * 1.0/2.0 * 1.0/r^6
            #        atoms.r12[indices[1],indices[2]] +=  4. * 1.0/2.0 * 1.0/r^12
            #    else 
            #        atoms.r6[indices[1],indices[2]] += 4. * 1.0/r^6
            #        atoms.r12[indices[1],indices[2]] += 4. * 1.0/r^12
            #    end
            end
        end
    end
#    println("Number of found interactions: ", foundInteractions)
    #return distMat
end

function set_random_unit_velocities!(atoms,nTypes, KE_max)
    nAtoms = length(atoms.positions)
#    masses = vcat(atoms.masses...)  # Flatten the masses vector.
    # Get a set of random unit vectors in 3D space.
    unit_vecs = [SVector{3,Float64}(randn(3)) for i=1:nAtoms]
    for i=1:nAtoms
        unit_vecs[i] /= norm(unit_vecs[i])
    end
    # In 3D the magnitude of the velocity should follow a r^1/3N distribution
    # where N is the number of atoms.  Generate a magnitude from this distribution.
    mag = rand()^(1/(3 * nAtoms))
    velocities = similar(unit_vecs)
    for iAtom=1:nAtoms
        scaledMag = mag * sqrt(2 * KE_max/nAtoms / atoms.masses[iAtom])
        velocities[iAtom] = scaledMag * unit_vecs[iAtom]
        println("check norm: ",norm(velocities[iAtom]))
        println("mag: ", scaledMag)
    end
    atoms.velocities .= velocities
#    scaledMag = mag * sqrt.([2 * KE_max ./ x for x in atoms.masses])
#    velocities = [scaledMag[x][y] .* unit_vecs[x][y] for x in 1:length(types) for y in 1:length(x)]# * unit_vecs
#    return velocities
end


function set_masses(atoms::atoms,masses::Vector{Vector{Float64}})
    for i=1:length(atoms.positions)
        atoms.masses[i] = masses[i]
    end
end

function set_masses(atoms::atoms,masses::Float64)
    for i=1:length(atoms.positions)
        atoms.masses[i] = masses
    end
end



function set_velocities(atoms::atoms,velocities::Float64)
    for i=1:length(atoms.positions)
        atoms.velocities[i] = SVector{3,Float64}(zeros(3))
    end
end
function precalc_LJ_distances!(atoms::atoms,cutoff)
    CartesianToDirect!(atoms)
    
    eVals = eigvals(transpose(atoms.latpar .* atoms.lVecs) * (atoms.latpar .* atoms.lVecs))
    maxN = cutoff/sqrt(minimum(eVals))
    
    loopBounds = SVector{3,Int64}(convert.(Int64,cld.(cutoff ,SVector{3,Float64}(norm(x) for x in eachcol(atoms.latpar * atoms.lVecs)) )))
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtom) in enumerate(atoms.positions)
        centerAtomType = atoms.atomTypes[iCenter]
        centerAtomC = DirectToCartesian(atoms.latpar * atoms.lVecs,centerAtom)
        get_single_atom_r6_r12!(atoms,cutoff,centerAtomC,centerAtomType,loopBounds)    # Find the contribution to the LJ energy for this centering atom.
        
    end
end

function eval_PE_LJ(atoms::atoms,model;force_recalc = false)
    if all(atoms.lj_vec .== 0.0)
 #       println("Calculating Distances")
        precalc_LJ_distances!(atoms,model.cutoff)
    elseif force_recalc
        atoms.lj_vec .= 0.0
        precalc_LJ_distances!(atoms,model.cutoff)
    end

#    println("vecs",atoms.lj_vec)
    coeffs = vcat(-1.0 * model.ϵ .* model.σ.^6, model.ϵ .* model.σ.^12)
    #println(coeffs)
    #println(atoms.lj_vec)
    totalEnergy = atoms.lj_vec' * coeffs

    return (totalEnergy + model.offset) * model.stdEnergy + model.meanEnergy
end

    


function gradientForce(model,atoms,atom,loopBounds; eps = 1e-3)
    fVec = zeros(3)
    # Find x-component of force
#    println("before")
#    display(atoms.positions[atom[1]][atom[2]])
    DirectToCartesian!(atoms)
    atoms.positions[atom] -= SVector{3,Float64}(eps,0,0)
    energyOne = eval_energy(atoms,model,force_recalc = true)
#    println("energyOne", energyOne)
    DirectToCartesian!(atoms)
    atoms.positions[atom] += SVector{3,Float64}(2 * eps,0,0)
    energyTwo = eval_energy(atoms,model,force_recalc = true)
#    println("energyTwo", energyTwo)
    DirectToCartesian!(atoms)
    atoms.positions[atom] -= SVector{3,Float64}(eps,0,0)  # Put it back where it was.
    fVec[1] = -(energyTwo - energyOne)/(2 * eps)
#    println("energies after 1")
#    display(energyTwo)
#    display(energyOne)
    # Find y-component of force
    atoms.positions[atom] -= SVector{3,Float64}(0,eps,0)
    energyOne = eval_energy(atoms,model,force_recalc = true)
#    println("energyOne", energyOne)
    DirectToCartesian!(atoms)
    atoms.positions[atom] += SVector{3,Float64}(0,2 * eps,0)
    energyTwo = eval_energy(atoms,model,force_recalc = true)
#    println("energyTwo", energyTwo)
    DirectToCartesian!(atoms)
    atoms.positions[atom] -= SVector{3,Float64}(0,eps,0)  # Put it back where it was.
    fVec[2] = -(energyTwo - energyOne)/(2 * eps)

 #   println("energies after 2")
 #   display(energyTwo)
 #   display(energyOne)
    # Find z-component of force
    atoms.positions[atom] -= SVector{3,Float64}(0,0,eps)  # Move atom downward
    energyOne = eval_energy(atoms,model,force_recalc = true)        # Calculate energy
#    println("energyOne", energyOne)
    DirectToCartesian!(atoms)
    atoms.positions[atom] += SVector{3,Float64}(0,0,2 * eps)  # Move atom upward
    energyTwo = eval_energy(atoms,model,force_recalc = true)           # Calculate energy
#    println("energyTwo", energyTwo)
    DirectToCartesian!(atoms)
    atoms.positions[atom] -= SVector{3,Float64}(0,0,eps)  # Put it back where it was.
    fVec[3] = -(energyTwo - energyOne)/(2 * eps)      # Calculate component of gradient
    DirectToCartesian!(atoms)
 #   println("energies after 3")
 #   display(energyTwo)
 #   display(energyOne)

    return fVec
end


function get_random_displacements(n_vecs,step_size)
    displacements = [zeros(SVector{3}) for i = 1:n_vecs]
    for j = 1:length(n_vecs)
        randVel = convert(SVector{3},(  2 * rand(3) .- 1.0))
        displacements[j] = step_size*randVel/norm(randVel)
    end
    return displacements


end


function eval_forces(atoms,model)
    forces = zeros(SVector{3,Float64},length(atoms.positions))

    if lowercase(model.name) == "lj"
        for (iAtom,atom) in enumerate(atoms.positions)  #Loop over the different atom types.
            forces[iAtom] = get_single_atom_force_lj(model,atoms,iAtom,@SVector[2,2,2])
        end

    end 
    return forces
end

function eval_forces!(atoms,model,forces)
#    forces = zeros(SVector{3,Float64},length(atoms.positions))

    if lowercase(model.name) == "lj"
        for (iAtom,atom) in enumerate(atoms.positions)  #Loop over the different atom types.
            forces[iAtom] = get_single_atom_force_lj(model,atoms,iAtom,@SVector[2,2,2])
        end

    end 
end

# Galilean Monte Carlo
function do_GMC!(config::atoms,model,walk_params,E_max)
    initialConfig = deepcopy(config)
    oldEnergy = ase.eval_energy(config,model)
    ase.DirectToCartesian!(config)
   # println("energy at start")
   # println(oldEnergy)
   # println(config.coordSys)

    displacements = get_random_displacements(length(config.positions),walk_params.MC_atom_step_size)
    forces = [zeros(SVector{3}) for x = 1:length(config.positions)]
    #display(velocities)
    n_reflect = 0
    n_reverse = 0
    n_accept = 0
    begin_energy = eval_energy(config,model,force_recalc = true)
    #println("Energy at beginning:", begin_energy)
    for iWalk = 1:walk_params.atom_traj_length
        energy_at_start_of_loop = eval_energy(config,model,force_recalc = true)  # Calculate the new energies
       # println("Energy at start of loop: ", energy_at_start_of_loop)
        last_good_positions = deepcopy(config.positions)
        last_good_displacements = deepcopy(displacements)
        config.positions .+= displacements  # Propogate the postions forward
        ase.mapIntoCell!(config)
        if any([isnan(y) for x in config.positions for y in x])
            println(displacements)
            error("Found NaNs")
        end
        energy_after_displacement = eval_energy(config,model,force_recalc = true)  # Calculate the new energies
        #println("Energy after displacement: ", energy_after_displacement)
        ase.DirectToCartesian!(config)
        if energy_after_displacement > E_max  # If we went uphill, we need to try and re-direct the displacements in the direction of the net force.
            #n_reflect += 1
            eval_forces!(config,model,forces)
#            for (iAtom,atom) in enumerate(config.positions)  #Loop over the different atom types.
#                forces[iAtom] = ase.get_single_atom_force_lj(model,config,iAtom,@SVector[2,2,2])
#            end
            # Sometimes the force vectors are zero, and therefore the normalization fails (divide by zero)
            # If the magnitude of the force is zero, then we don't need to do anything.
            nHat = [norm(x) > 0.0 ? x / norm(x) : x for x in forces]
            #nHat = [x ./ norm.(x) for x in forces]
            for i = 1:length(config.positions)
                @inbounds displacements[i] -= (2 * (displacements[i]' * nHat[i])) * nHat[i]
            end
            if any([isnan(z) for x in displacements for y in x for z in y])
                println(displacements)
                println(nHat)
                println(forces)
                error("Found NaNs in displacement 1")
            end
            config.positions += displacements
            n_reflect += 1
            ase.mapIntoCell!(config)
            energy_after_deflection = ase.eval_energy(config,model,force_recalc = true)
         #   println("Energy after deflection: ", energy_after_deflection)
            if energy_after_deflection > E_max  # Reverse
          #      println("Reversing")
                config.positions = deepcopy(last_good_positions)
                displacements = -1.0 * deepcopy(last_good_displacements)
                if any([isnan(z) for x in displacements for y in x for z in y])
                    println(displacements)
                    error("Found NaNs in displacement 2")
                end
                n_reflect -= 1
                n_reverse += 1
                # If the very last iteration was a failed reverse, then we don't have the total energy stored anywhere.
            else
           #     println("deflection accepted--------------------------------------------------------")
            #    println("Energy:", eval_energy(config,model,force_recalc = true))
                config.model_energy = energy_after_deflection
            end
        else
           # println("Displacement accepted")
           # println("Energy:", eval_energy(config,model,force_recalc = true))
            n_accept += 1
            config.model_energy = energy_after_displacement
        end
        

    end

#    println("No of reflections: " ,n_reflect)
#    println("No of reverses: ",n_reverse)
#    println("No of straight acceptances: " ,n_accept)
#    println("Beginning Energy: ", begin_energy)
#    println("End energy: ", config.model_energy )
#    println("Energy Threshold: ", E_max)
    return n_reflect + n_reverse, n_reflect#    if newEnergy > oldEnergy
#        return initialConfig
#    else
#        return config
#    end
end

function get_nTypes(atoms)
return [count(==(atomT),atoms.atomTypes) for atomT in unique(atoms.atomTypes)]
end

# Perform a random walk on a single configuration subject to the constraint that the total energy be less than the cutoff
# Not walking using the force vector (GMC), just random walks.  This may not work well for some systems.
function do_MC!(atoms::atoms,model,walk_params,E_max) 
    # Loop over the number of random walks to take.
    n_accept = 0
    for iWalk in 1:walk_params.atom_traj_length
        #@printf "Step in random walk %3i\n" iWalk
        #Loop over all of the atoms in the simulation.
        for (iAtom,atom) in enumerate(atoms.positions)
            #@printf "Atom being moved. Type: %3i Number: %3i\n" iType iAtom 
            # Get a random displacement vector
            randDisplace = (2*rand(3) .- 1.0)*walk_params.MC_atom_step_size
#            println("before")
#            display(atoms.positions[iType][iAtom])
            atoms.positions[iAtom] .+= randDisplace
            mapIntoCell!(atoms)
#            println("after")
#            display(atoms.positions[iType][iAtom])
#            println(atoms.coordSys)
            newEnergy = eval_energy(atoms,model) # Want total energies for NS
            # If the move resulted in a higher energy, undo the move and go to the next atom.
            if newEnergy > E_max  # reject move
                #println("Rejected")
                atoms.positions[iAtom] .-= randDisplace
            else # Otherwise, update the energy 
                #println("Accepted")
                atoms.model_energy = newEnergy
                n_accept += 1
            end     
        end   
#        randDisplacement = [[ (rand(3).-0.5)*0.1 for i =1:atoms.nType[j]] for j = 1:atoms.order]
#        atoms.positions .+= randDisplacement  # Use the displacement to move this atom.

       # println(totalEnergy(atoms,model))
    end
    return (n_steps,n_accept)
end

function get_momenta(atoms)
    return atoms.masses .* atoms.velocities

end

function eval_KE(atoms)
    KE = 1/2 * atoms.masses .* [x' * x for x in atoms.velocities]
    return sum(KE)

end


function eval_energy(atoms,model;P = 0,do_KE = true, force_recalc = false)
#    println("here", P * cell_volume(atoms))
#    println("here", eval_KE(atoms))
    energy = P * cell_volume(atoms)
    if do_KE
        energy += eval_KE(atoms)
    end
#    println("here", energy)
 #   println(typeof(model))
    if lowercase(model.name) == "lj"
        energy +=  eval_PE_LJ(atoms,model,force_recalc = force_recalc)
    else
        error("I don't know what model to use.")
    end
#    println(energy)
#    if typeof(model) == LennardJones.model{Matrix{Float64}}
#        println("Found Lennard Jones model")
#        return LennardJones.totalEnergy(atoms,model)
#
#
#        # Add other use cases as models become available
#
#    end
    return energy
end

function minkowski_reduce(ABC,eps)
    #ABC = hcat(A,B,C)
    outABC = deepcopy(ABC)
    limit = 10
    for it=1:limit
        norms = [norm(x) for x in eachcol(outABC)]
        outABC = outABC[:,sortperm(norms)]
        outABC = reduce_C_in_ABC(outABC,1e-7)
        if norm(outABC[:,3]) >= norm(outABC[:,2]) - eps
            break
        end
    end

    if !minkowski_check(outABC,eps)
        error("Minkowski conditions not met")
    end
    if det(outABC') < 0
        save = outABC[:,2]
        outABC[:,2] .= outABC[:,3]
        outABC[:,3] .= save
    end
    return outABC
end


function minkowski_check(ABC,eps)

    b1 = ABC[:,1]
    b2 = ABC[:,2]
    b3 = ABC[:,3]

    minkowski_check = true

    if norm(b1) > norm(b2) + eps
        minkowski_check = false
        println("Minkowski_condition 1 failed: b1 > b2")
    end

    if norm(b2) > norm(b3) + eps
        minkowski_check = false
        println("Minkowski_condition 2 failed: b2 > b3")
    end

    if norm(b2) > norm(b1 + b2) + eps
        minkowski_check = false
        println("Minkowski_condition 3 failed: b2 > b1 + b2")
    end

    if norm(b2) > norm(b1-b2) + eps
        minkowski_check = false
        println("Minkowski_condition 4 failed: b2 > b1-b2")
    end

    if norm(b3) > norm(b1+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 5 failed: b3 > b1+b3")
    end

    if norm(b3) > norm(b3-b1) + eps
        minkowski_check = false
        println("Minkowski_condition 6 failed: b3 > b3-b1")
    end

    if norm(b3) > norm(b2+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 7 failed: b3 > b2+b3")
    end

    if norm(b3) > norm(b3-b2) + eps
        minkowski_check = false
        println(b3)
        println(b2)
        println("Minkowski_condition 8 failed: b3 > b3-b2")
    end

    if norm(b3) > norm(b1+b2+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 9 failed: b3 > b1+b2+b3")
    end

    if norm(b3) > norm(b1-b2+b3) + eps
        minkowski_check = false
        println("Minkowski_condition 10 failed: b3 > b1-b2+b3")
    end

    if norm(b3) > norm(b1+b2-b3) + eps
        minkowski_check = false
        println("Minkowski_condition 11 failed: b3 > b1+b2-b3")
    end

    if norm(b3) > norm(b1-b2-b3) + eps
        minkowski_check = false
        println("Minkowski_condition 12 failed: b3 > b1-b2-b3")
    end

    return minkowski_check
end






function vegardsVolume(elements,atom_counts,volume)
    if nothing in [findfirst(x-> x == lowercase(i),lowercase.(keys(element_volume))) for i in elements]
        return 1.0
    end
    nAtoms = sum(atom_counts)
    nTypes = length(atom_counts)
    concentrations = [x/nAtoms for x in atom_counts]
    if length(elements) != nTypes
        error("Not enough elements specified")
    end

    #  The next two lines are Wiley's way.
    # Find lat pars if crystal was pure.
    purelatPars = [( element_volume[elements[x]]/(volume/nAtoms) )^(1. /3.) for x =1:nTypes]
    # Concentration-weighted average over lattice parameters
    wiley = sum([purelatPars[x] * concentrations[x] for x =1:nTypes])


    # The next two lines are my way.
    # concentration-weighted average of volumes
    avgVolume = sum([element_volume[elements[x]] *concentrations[x] for x =1:nTypes] )
    # Calculate lattice parameter
    mine = ( avgVolume/(volume/nAtoms) )^(1/3.)
#    print("My approach: {}.  Wiley's approach: {}".format(mine,wiley))
    return wiley
end


function do_MD!(atoms::atoms, model, walk_params;E_cutoff = NaN,KE_cutoff = NaN)
    # Just copied over from MD. Needs checked!!!

    N = length(atoms.positions)
    forces = zeros(SVector{3,Float64},N)
    pre_positions = deepcopy(atoms.positions)
    pre_velocities = deepcopy(atoms.velocities)
    idx = 1
    begin_energy = eval_energy(atoms,model,force_recalc = true)
    begin_KE = eval_KE(atoms)
    for i=1:walk_params.atom_traj_length
        eval_forces!(atoms,model,forces)
#        println("check energy: ", eval_energy(atoms,model,force_recalc = true))
#        for n=1:N
#            forces[n] = get_single_atom_force_lj(model,atoms,n,[2,2,2])
#        end
#        display(forces)
        atoms.velocities .+= forces ./ atoms.masses * 0.5 * walk_params.MD_time_step
 #       println("Check velocities")
  #      display(atoms.velocities)
        atoms.positions .+= atoms.velocities * walk_params.MD_time_step
        mapIntoCell!(atoms)  # We may have gone outside the cell, put the atoms back
        eval_forces!(atoms,model,forces)

##        for n=1:N
#            forces[n] = get_single_atom_force_lj(model,atoms,n,[2,2,2])
#        end
        atoms.velocities .+= forces ./ atoms.masses * 0.5 * walk_params.MD_time_step
  #      display(forces)
    end 
    final_energy = eval_energy(atoms,model, force_recalc = true)
    final_KE = eval_KE(atoms)
   # println("Begin Energy:-------------------------------------> ", begin_energy)
   #println("Final Energy: ", final_energy)
    if abs(final_energy - begin_energy) > walk_params.MD_reject_eps * final_KE
#        println("Rejecting on non-energy conserving MD walk")
        n_accept = 0
        atoms.positions .= deepcopy(pre_positions)
        atoms.velocities .= deepcopy(pre_velocities)
 #       println("Beginning Energy:", begin_energy)
 #       println("Final Energy:", final_energy)
        return (1,n_accept)
    elseif !isnan(E_cutoff) && !isnan(KE_cutoff)
        if final_energy > E_cutoff || final_KE > KE_cutoff
  #          println("rejecting")
            #println("End energy (check): ", final_energy)
            #println("E-max (check): ", E_max)
            n_accept = 0
            atoms.positions .= deepcopy(pre_positions)
            atoms.velocities .= deepcopy(pre_velocities)
        else#if final_energy < E_cutoff || final_KE < KE_cutoff
   #         println("Accepting")
            n_accept = 1
            atoms.model_energy = final_energy
            atoms.velocities *= -1.0
        end

        return (1,n_accept)
    else
        println(E_cutoff, KE_cutoff)
        error("I don't know what I'm suppose to do. I thought one of the other two branches would catch" )
        
    end
end


end