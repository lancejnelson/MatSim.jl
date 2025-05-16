module nSampling

# Import local libraries
using ase
using LennardJones

# Import packaged libraries
using DataSets:DataSet
using StaticArrays
using LinearAlgebra
using Distributions
using Random


struct NS_walker_params
    n_single_walker_steps:: Int
    n_atom_steps:: Int
    n_cell_volume_steps:: Int
    n_cell_shear_steps:: Int
    n_cell_stretch_steps:: Int
    atom_traj_length:: Int
    min_aspect_ratio:: Float64
    max_volume_per_atom:: Float64
    volume_step_size:: Float64
    shear_step_size:: Float64
    stretch_step_size:: Float64
    atom_algorithm:: String
end


struct NS
    n_walkers:: Int
    n_cull:: Int
    n_iter:: Int
    cell_P:: Float64
    eps:: Float64
    walker_params:: NS_walker_params
    walkers::Vector{ase.atoms}
end



# Initializes a Nested Sampling run, taking in the options, atom types, and model parameters. (Currently set up for Lennard Jones)
function initialize(inputs,species,model)
    # Add all atoms of all species types. inputs is a list of integers corresponding to each species.
    nAtoms = sum(inputs["n_Atoms"])
    
    # Get a random lattice constant spacing depending on the number of atoms and the maximum volume allowed each atom
    latticeConstant = (nAtoms * inputs["max_volume_per_atom"] * rand()^(1/(nAtoms + 1)))^(1/3)
    lVecs = SMatrix{3,3,Float64,9}(latticeConstant * [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

    # Initialize a lattice for each walker
    configs = [ase.initialize_cell_shape(latticeConstant,lVecs) for i in 1:inputs["n_walkers"]]
#    parsedlVecs = [parse.(Float64,split(x[1])) for x in inputs["lVecs"]]
#    lVecs = hcat(parsedlVecs...)
#    println("lattice constant: ", latticeConstant)
#    println("lattice vectors: ", lVecs)
#    configs = [ase.initialize_with_random_positions(latticeConstant,lVecs,inputs["n_Atoms"],inputs["minSep"],species) for i in 1:inputs["n_walkers"]]
    
    # Load in the walker parameters for all walkers
    # TODO Should this be added to the NS.YAML or something similar?
    walk_params = Dict("n_steps" => 10,
                       "shear_step_size" => inputs["shear_step_size"],
                       "stretch_step_size" => inputs["stretch_step_size"])
    # Now do a random walk on the cell shape and place atoms at random positions
    for i in configs
        #println("walker: ", i)
        do_cell_shape_walk!(i,walk_params)
        ase.set_atoms_random!(i,inputs["n_Atoms"],inputs["minSep"],species)
        ase.set_masses(i,1.0)
        ase.set_random_unit_velocities!(i,nAtoms,2.0)
        i.model_energy = ase.eval_energy(i,model)
    end
    walker_params = NS_walker_params(inputs["n_single_walker_steps"],inputs["n_atom_steps"],
        inputs["n_cell_volume_steps"],inputs["n_cell_shear_steps"],inputs["n_cell_stretch_steps"],
        inputs["atom_traj_length"],inputs["min_aspect_ratio"],inputs["max_volume_per_atom"],
        inputs["volume_step_size"],inputs["shear_step_size"],inputs["stretch_step_size"],inputs["atom_algorithm"])
    return NS(inputs["n_walkers"],inputs["n_cull"],inputs["n_iter"],inputs["cell_P"],inputs["eps"],
              walker_params,configs)
end


function run_NS(NS::NS,LJ::LennardJones.model)

    V = (NS.n_walkers - NS.n_cull + 1)/(NS.n_walkers + 1)
    cDir = pwd()

    i = 1
    io = open(joinpath(cDir,"NS.out"),"w")
    write(io, "N_walkers = " * string(NS.n_walkers) * "\n")
    write(io, "N_cull = " * string(NS.n_cull) * "\n")
    write(io, "N_steps_per_walker = " * string(NS.n_iter) * "\n")
    write(io, "eps = " * string(NS.eps) * "\n")
    
    while V > NS.eps
        println(i)
        println("V = ", V)
        ## Find the top Kr highest energies
        sEnergies =  reverse(sortperm([i.energyPerAtomModel for i in NS.walkers]))
       # println(sort([i.energyPerAtomModel for i in NS.walkers]))
        energyCutoff = NS.walkers[sEnergies[NS.n_cull]].energyPerAtomModel  # Find the 
        # Which configs need to be thrown out.
        forDelete = sEnergies[1:NS.n_cull]
        # Which configs can be kept
        keeps = sEnergies[NS.n_cull + 1: end]

        println("Energy cutoff")
        display(energyCutoff)
        for replace_walker in forDelete
            #Copy one of the configs that didn't get thrown out as the starting point
            #println("Initializing random walk to replace configuration ", i)
            #display(sample(keeps))
            NS.walkers[replace_walker] = NS.walkers[sample(keeps,1)[1]]
            walk_single_walker!(NS.walkers[replace_walker],LJ,NS.walk_params)
#            randomWalk!(NS.walkers[replace_walker],LJ,energyCutoff,NS.n_iter)
        end
        i += 1
        V = ((NS.n_walkers - NS.n_cull + 1)/(NS.n_walkers + 1))^i
        write(io,string(V) * " ")
        write(io,string(energyCutoff) * " \n")

    end

    close(io)
end


function walk_single_walker!(atoms::ase.atoms, model, walk_params::NS_walker_params,E_max)
    possible = [do_cell_volume_step, do_cell_shear_step, do_cell_stretch_step, do_atoms_step]
    # Loop over the number of random walks to take.
    for iWalk in 1:walk_params.n_single_walker_steps
#        println("Current Energy", ase.eval_energy(atoms,model))
        move = rand(1:length(possible))
        possible[move](atoms,model,walk_params,E_max)

    end

end

function do_atoms_step(atoms::ase.atoms, model,walk_params,E_max) #atoms::ase.atoms, model,walk_params,E_max)
    if walk_params.atom_algorithm == "GMC"
        ase.do_GMC!(atoms,walk_params.n_atom_steps,model,E_max)
    elseif walk_params.atom_algorithm == "MC"
        ase.do_MC!(atoms,walk_params.n_atom_steps,model,E_max)
    elseif walk_params.atom_algorithm == "MD"
        ase.do_MD!(atoms,walk_params.n_atom_steps,model,E_max)
    end
end

function do_cell_shear_step(atoms::ase.atoms, model,walk_params,E_max)
    (p,T) = propose_shear_step(atoms,walk_params.shear_step_size)
    do_cell_step!(atoms,p,T,E_max,walk_params,model)

end

function do_cell_stretch_step(atoms::ase.atoms, model,walk_params,E_max)
    (p,T) = propose_stretch_step(atoms,walk_params.stretch_step_size)
    do_cell_step!(atoms,p,T,E_max,walk_params,model)

end


function do_cell_volume_step(atoms::ase.atoms, model,walk_params,E_max)
    (p,T) = propose_volume_step(atoms,walk_params.volume_step_size)
    do_cell_step!(atoms,p,T,E_max,walk_params,model)

end

function do_cell_step!(atoms,p_accept,T,E_max,ns_params,model)
    println("Doing cell step")
    if (p_accept < 1.0) & (rand() > p_accept)
        println("rejected")
        return false   
    end

    newCell = atoms.lVecs * T
    newVolume = ase.cell_volume(newCell)
    if newVolume > ns_params.max_volume_per_atom * atoms.nAtoms
 #       println("Volume too large")
 #       println("new volume: ", newVolume)
 #       println("old volume: ", ase.cell_volume(atoms.lVecs))
        return false
    end
    if ase.min_aspect_ratio(newCell) < ns_params.min_aspect_ratio

#        println("Aspect ratio too small")
#        println("new aspect ratio: ", ase.min_aspect_ratio(newCell))
#        println("old aspect ratio: ", ase.min_aspect_ratio(atoms.lVecs))
        return false
    end

    oldCell = atoms.lVecs
    oldPositions = atoms.positions
    println("before")
    display(ase.eval_energy(atoms,model))
    
    atoms.lVecs = newCell
    ase.set_cell!(atoms,newCell,scale_atoms = true)
    println("after ")
    display(ase.eval_energy(atoms,model))
    newEnergy = ase.eval_energy(atoms,model)
    if newEnergy < E_max
        println("Accepted")
        atoms.fitEnergy = newEnergy
    else
        ase.set_cell!(atoms,oldCell)
        atoms.positions = oldPositions # Can I just set cell with scale_atoms = true again?
    end

end


function propose_volume_step(atoms::ase.atoms, stepsize::Float64)
    dV = rand(Normal(0,stepsize * atoms.nAtoms))
    # Get the current volume

    orig_V = ase.cell_volume(atoms)
    # Calculate the new volume
    new_V = orig_V + dV
    if new_V < 0
        println("Negative volume encountered, resetting to positive value.")
        new_V = abs(new_V)
    end
    # Calculate the scale factor
    transform = UniformScaling((new_V / orig_V)^(1.0/3.0))
    
    # This skews the probability towards larger volumes..
    # Because if new_V > orig_V, then the probability of accepting the move is 1.0, if new_V < orig_V, then the probability of accepting the move is less than 1.0. (Note... I'm still not sure why we want to do this. Waiting on a response from Livia)
    p_accept = minimum([1.0, (new_V/orig_V)^atoms.nAtoms])
    
    
    return (p_accept, transform)
end


function propose_shear_step(atoms::ase.atoms,stepsize::Float64)
    rmv_vec = rand(1:3)
    remaining_vecs = [x for x in 1:3 if x != rmv_vec]
    new_cell = convert(Matrix,copy(atoms.lVecs))
    v1 = atoms.lVecs[:,remaining_vecs[1]]
    v2 = atoms.lVecs[:,remaining_vecs[2]]
    v1 /= sqrt(v1' * v1)
    v2 -= v1 * (v1' * v2)  # Make v2 orthogonal to v1
    v2 /= sqrt(v2' * v2)
    if abs(v1' * v2) > 1e-4
        println("v1 and v2 are not orthogonal, something is wrong.")
    end
    rv1 = rand(Normal(0,stepsize))
    rv2 = rand(Normal(0,stepsize))
    new_cell[:,rmv_vec] = atoms.lVecs[:,rmv_vec] + rv1 * v1 + rv2 * v2
    transform = new_cell * inv(atoms.lVecs)
    check = transform * atoms.lVecs
    return (1.0, transform)  # Return the shear matrix and the probability of accepting the move.
end


function propose_stretch_step(atoms::ase.atoms,stepsize::Float64) 
    # Yes, I know that I don't need an atoms object to do this, but to stay consistent with the shear step, I'm going to keep it.


#    a = stepsize .* (2 .* rand(2) .- 1)  # Get 3 random floats
#    A = [a[1] 0  0
#         0  a[2] 0
#         0  0  1/a[1]/a[2]]
#    return A

    rmv_vec_one = rand(1:3)
    rmv_vec_two = rand(1:3)
    if rmv_vec_one == rmv_vec_two
        #println("found duplicate vectors, regenerating rmv_vec_two")
        rmv_vec_two = (rmv_vec_two + 1) % 3 + 1
    end
    remaining_vec = [x for x in 1:3 if x != rmv_vec_one && x != rmv_vec_two][1]

    rv = rand(Normal(0,stepsize))
    transform = zeros(3,3)
    transform[remaining_vec,remaining_vec] = 1.0
    transform[rmv_vec_one,rmv_vec_one] = exp(rv)
    transform[rmv_vec_two,rmv_vec_two] = exp(-rv)
    return (1.0, transform)
end



function do_cell_shape_walk!(atoms::ase.atoms, walk_params::Dict)
    possibilities =[propose_shear_step,propose_stretch_step]

    for i =1:walk_params["n_steps"]
        shuffle!(possibilities)
        for j in possibilities
            p_accept, transform = j(atoms,walk_params[string(j)[9:end] * "_size"])
            if rand() < p_accept
                atoms.lVecs = transform * atoms.lVecs
                ase.DirectToCartesian!(atoms)
           #     println("Accepted")
            end
        end
    end

end

function get_random_displacements(n_vecs)
    displacements = [zeros(SVector{3}) for i = 1:n_vecs]
    for j = 1:length(n_vecs)
        randVel = convert(SVector{3},(  2 * rand(3) .- 1.0))
        displacements[j] = 0.2*randVel/norm(randVel)
    end
    return displacements


end






end