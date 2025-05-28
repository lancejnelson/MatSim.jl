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


mutable struct NS_walker_params
    n_single_walker_steps:: Int
    MC_atom_step_size:: Float64
    MD_time_step:: Float64
    MD_reject_eps:: Float64
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
    KE_max:: Float64
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

function NS_walker_params(inputs)
    final_settings = Dict()
    for (key,value) in NS_defaults()
        if !(key in keys(inputs))
            # Use the default
            final_settings[key] = value
        else
            # Use custom setting
            final_settings[key] = inputs[key]
        end

    end
    nAtoms = sum(inputs["n_Atoms"])
    max_lc = (nAtoms * inputs["max_volume_per_atom"])^(1/3)  # This is the maximum lattice constant

     return NS_walker_params(final_settings["n_single_walker_steps"],max_lc * final_settings["MC_atom_step_size"],final_settings["MD_time_step"],final_settings["MD_reject_eps"],final_settings["n_atom_steps"],
        final_settings["n_cell_volume_steps"],final_settings["n_cell_shear_steps"],final_settings["n_cell_stretch_steps"],
        final_settings["atom_traj_length"],final_settings["min_aspect_ratio"],final_settings["max_volume_per_atom"],
        final_settings["volume_step_size"],final_settings["shear_step_size"],final_settings["stretch_step_size"],inputs["KE_max"],final_settings["atom_algorithm"])
end

function NS_defaults()
    defaults = Dict("min_atom_separation"=> 1.0,
                    "n_walkers" => 100, 
                    "n_cull" => 1, 
                    "atom_traj_length"=> 8, 
                    "eps"=> 1e-6,
                    "n_iter" => 2000,
                    "atom_algorithm"=> "MD", 
                    "MC_atom_step_size" => 0.5, 
                    "MD_time_step"=>0.1,
                    "MD_reject_eps"=>1e-2,
                    "n_single_walker_steps"=> 600, 
                    "n_atom_steps"=> 100, 
                    "n_cell_volume_steps"=> 5, 
                    "n_cell_shear_steps"=>8,
                    "n_cell_stretch_steps"=>8, 
                    "cell_P"=>0.0, 
                    "min_aspect_ratio" => 0.8, 
                    "max_volume_per_atom"=>32.0,
                    "volume_step_size" => 0.5, 
                    "shear_step_size"=>0.1, 
                    "stretch_step_size"=>0.1)

    return defaults
end



function initialize(inputs,model)
    nAtoms = sum(inputs["n_Atoms"])
    max_lc = (nAtoms * inputs["max_volume_per_atom"])^(1/3)  # This is the maximum atom separation.
    latticeConstant = (nAtoms * inputs["max_volume_per_atom"] * rand()^(1/(nAtoms + 1)))^(1/3)
    lVecs = SMatrix{3,3,Float64,9}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    
#    error("stopping")
    configs = [ase.initialize_cell_shape(latticeConstant,lVecs,nAtoms) for i in 1:inputs["n_walkers"]]
    if haskey(inputs, "cell_P")
        inputs["KE_max"] = 3/2 * inputs["cell_P"] * inputs["max_volume_per_atom"] * nAtoms
    end
    walker_params = NS_walker_params(inputs)
#    walk_params = Dict("n_steps" => 10,
#                       "shear_step_size" => inputs["shear_step_size"],
#                       "stretch_step_size" => inputs["stretch_step_size"])
    # Now do a random walk on the cell shape and place atoms at random positions
    save_n_steps = walker_params.n_single_walker_steps
    walker_params.n_single_walker_steps = 10
    for i in configs
        do_cell_shape_walk!(i,walker_params)
        ase.set_atoms_random!(i,inputs["n_Atoms"],inputs["min_atom_separation"],10 * model.cutoff,model.system)
        ase.set_masses(i,1.0)
        if lowercase(inputs["atom_algorithm"]) == "md"
            ase.set_random_unit_velocities!(i,nAtoms,walker_params.KE_max)
        else
            ase.set_velocities(i,0.0)
        end
        i.model_energy = ase.eval_energy(i,model,force_recalc = true)
    end
    walker_params.n_single_walker_steps = save_n_steps
    return NS(inputs["n_walkers"],inputs["n_cull"],inputs["n_iter"],inputs["cell_P"],inputs["eps"],
              walker_params,configs)
end

function adjust_step_sizes!(walk_params,key,rate)

    if rate > 1//2
        if key == "volume"
            walk_params.volume_step_size *=1.05
        elseif key == "shear"
            walk_params.shear_step_size *=1.05
        elseif key == "stretch"
            walk_params.stretch_step_size *=1.05
        elseif key == "atoms"
            if walk_params.atom_algorithm == "MD"
                walk_params.MD_time_step *=1.05
            else
                walk_params.MC_atom_step_size *=1.05
            end
        end
        return false
    elseif rate < 1//4
        if key == "volume"
            walk_params.volume_step_size *=0.95
        elseif key == "shear"
            walk_params.shear_step_size *=0.95
        elseif key == "stretch"
            walk_params.stretch_step_size *=0.95
        elseif key == "atoms"
            if walk_params.atom_algorithm == "MD"
                walk_params.MD_time_step *=0.95
            else
                walk_params.MC_atom_step_size *=0.95
            end
        end
        return false
    else
        return true
    end

end

function tune_step_sizes!(NS,model::LennardJones.model)
    sEnergies =  reverse(sortperm([i.model_energy for i in NS.walkers]))
    keeps = sEnergies[NS.n_cull + 1: end]
    E_max = NS.walkers[sEnergies[1]].model_energy
    atoms = NS.walkers[sample(keeps,1)[1]]
   
    allGood = false
    index = 0
    orig_atoms = deepcopy(atoms)
    while !allGood
        index += 1
        a_rates = walk_single_walker!(atoms,model,NS.walker_params,E_max)
        display(a_rates)
        display(NS.walker_params)
        rates_good = []
        for key in keys(a_rates)
            rate = a_rates[key][2]/a_rates[key][1]
            push!(rates_good,adjust_step_sizes!(NS.walker_params,key,rate))

        end
        println("rates good")
        println(rates_good)
        allGood =  all(rates_good)
        println(allGood)
        atoms = deepcopy(orig_atoms)
        if index > 200
            error("Taking too long")
        end
    end
        

#    possible = [do_cell_volume_step, do_cell_shear_step, do_cell_stretch_step, do_atoms_step]
#    # Loop over the number of random walks to take.
#    n_try = walk_params.n_single_walker_steps
#    n_accept = 0
#    idx = 0
#    acceptance_rates = Dict("stretch"=>0//n_try,"shear"=>0//n_try,"volume"=>0//n_try, "atoms"=>0//n_try)
#
#    orig_config = deepcopy(atoms)
##    volume_accept_rate = 0//n_try
##    while volume_accept_rate > 1//2  || volume_accept_rate < 1//4
##        n_accept = 0
##        for i = 1:walk_params.n_single_walker_steps
##            (__,n_accept_iter) = do_cell_volume_step(atoms,model,walk_params,E_max )
##            n_accept += n_accept_iter
##        end
##        volume_accept_rate = n_accept/n_try
##        if volume_accept_rate > 1//2
##            walk_params.volume_step_size *= 1.05
##        elseif volume_accept_rate < 1//4
##            walk_params.volume_step_size *= 0.95
##        end
##        println("Current step size volume: ", walk_params.volume_step_size)
##        println("Current accept rate: ", volume_accept_rate)
##    end
##
##
##    shear_accept_rate = 0//n_try
##    while shear_accept_rate > 1//2  || shear_accept_rate < 1//4
##        n_accept = 0
##        for i = 1:walk_params.n_single_walker_steps
##            (__,n_accept_iter) = do_cell_shear_step(atoms,model,walk_params,E_max )
##            n_accept += n_accept_iter
##        end
##        shear_accept_rate = n_accept/n_try
##        if shear_accept_rate > 1//2
##            walk_params.shear_step_size *= 1.05
##        elseif shear_accept_rate < 1//4
##            walk_params.shear_step_size *= 0.95
##        end
##        println("Current step size shear: ", walk_params.shear_step_size)
##        println("Current accept rate: ", shear_accept_rate)
##    end
##
##    stretch_accept_rate = 0//n_try
##    while stretch_accept_rate > 1//2  || stretch_accept_rate < 1//4
##        n_accept = 0
##        for i = 1:walk_params.n_single_walker_steps
##            (__,n_accept_iter) = do_cell_stretch_step(atoms,model,walk_params,E_max )
##            n_accept += n_accept_iter
##        end
##        stretch_accept_rate = n_accept/n_try
##        if stretch_accept_rate > 1//2
##            walk_params.stretch_step_size *= 1.05
##        elseif stretch_accept_rate < 1//4
##            walk_params.stretch_step_size *= 0.95
##        end
##        println("Current step size stretch: ", walk_params.stretch_step_size)
##        println("Current accept rate: ", stretch_accept_rate)
##    end    
#
#
#
#    atom_accept_rate = 0//n_try
#    index = 0
#    while atom_accept_rate > 1//2  || atom_accept_rate < 1//4
#        index += 1
#        n_accept = 0
#        for i = 1:walk_params.n_single_walker_steps
#            (__,n_accept_iter) = do_atoms_step(atoms,model,walk_params,E_max, KE_max = walk_params.KE_max )
#            n_accept += n_accept_iter
#        end
#        atom_accept_rate = n_accept//n_try
#        if atom_accept_rate > 1//2
#            walk_params.MC_atom_step_size *= 1.05
#        elseif atom_accept_rate < 1//4
#            walk_params.MC_atom_step_size *= 0.95
#        end
#        atoms = deepcopy(orig_config)
#        println("Current step size atom: ", walk_params.MD_time_step)
#        println("Current accept rate: ", atom_accept_rate, "<--------------------------------------------------------->")
#        if index > 100
#            error("Too long")
#        end
#    end    
            #    for iWalk in 1:walk_params.n_single_walker_steps

   # for i =1:10 # How many steps should we take to do this?
   #     idx += 1
#  #      println("Current Energy", ase.eval_energy(atoms,model))
   #     move = rand(1:length(possible))
   #     (tried,accepted) = possible[move](atoms,model,walk_params,E_max)
   #     acceptance_rates[string(possible_move)[9:end-5]] += accepted//n_try
#  #      n_try += n_iter_try
#  #      n_accept += n_iter_accept
   #     println("Walk Type: ", move)
   #     println("Index: ", idx)
   #     println("n_accept: ", n_accept)
   #     println("n_try: ", n_try)
   #     if idx > 300
   #         error("Too long.. Stopping")
   #     end
   #     
   #     
   # end

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
    
#    while V > NS.eps
    for i= 1:100
        println(i)
        println("V = ", V)
        ## Find the top Kr highest energies
        sEnergies =  reverse(sortperm([i.model_energy for i in NS.walkers]))
       # println(sort([i.energyPerAtomModel for i in NS.walkers]))
        E_max = NS.walkers[sEnergies[NS.n_cull]].model_energy  # Find the 
        # Which configs need to be thrown out.
        forDelete = sEnergies[1:NS.n_cull]
        # Which configs can be kept
        keeps = sEnergies[NS.n_cull + 1: end]

        println("Energy cutoff")
        display(E_max)

        if i %25 == 0  # 25 is pretty arbitrary.. Need a better way to see if need to re-tune
            println("Stopping to retune step sizes")
            tune_step_sizes!(NS,LJ)
        end
#        display(E_max - ase.eval_KE(NS.walkers[sEnergies[NS.n_cull]]))
        for replace_walker in forDelete
            #Copy one of the configs that didn't get thrown out as the starting point
            #println("Initializing random walk to replace configuration ", i)
            #display(sample(keeps))
            NS.walkers[replace_walker] = deepcopy(NS.walkers[sample(keeps,1)[1]])
            walk_single_walker!(NS.walkers[replace_walker],LJ,NS.walker_params,E_max)
#            randomWalk!(NS.walkers[replace_walker],LJ,energyCutoff,NS.n_iter)
        end
        i += 1
        V = ((NS.n_walkers - NS.n_cull + 1)/(NS.n_walkers + 1))^i
        write(io,string(V) * " ")
        write(io,string(E_max) * " \n")

    end

    close(io)
end


function walk_single_walker!(atoms::ase.atoms, model, walk_params::NS_walker_params,E_max)
    possible = [do_cell_volume_step, do_cell_shear_step, do_cell_stretch_step, do_atoms_step]
    # Loop over the number of random walks to take.
    n_try = 0
    n_accept = 0
    idx = 0
    acceptance_rates = Dict("stretch"=>[0,0],"shear"=>[0,0],"volume"=>[0,0], "atoms"=>[0,0])

    for iWalk in 1:walk_params.n_single_walker_steps
#    while n_accept < walk_params.n_single_walker_steps
        idx += 1
#        println("Current Energy", ase.eval_energy(atoms,model))
        move = rand(1:length(possible))
        tried = 0
        moveType = split(string(possible[move]),"_")[end-1]
#        println("move type: -------------------------------------------->", moveType)
        before_energy = ase.eval_energy(atoms,model,force_recalc = true)
        (tried,accepted) = possible[move](atoms,model,walk_params,E_max)
        after_energy = ase.eval_energy(atoms,model,force_recalc = true)
 #       println("Energy before move: ", before_energy)
 #       println("Energy after move: ", after_energy)
 #       println("E-max: ", E_max)
 #       println("Tried: ", tried)
 #       println("Accepted: ", accepted)
        acceptance_rates[moveType][1] += tried
        acceptance_rates[moveType][2] += accepted
#        n_try += n_iter_try
#        n_accept += n_iter_accept
#        println("Walk Type: ", move)
#        println("Index: ", idx)
#        println("n_accept: ", n_accept)
#        println("n_try: ", n_try)
        if idx > 10000
            error("Too long.. Stopping")
        end
        
        
    end
#    for key in keys(acceptance_rates)
#        push!(acceptance_rates[key],acceptance_rates[key][2]/acceptance_rates[key][1])
#    end
    return acceptance_rates
end

function do_atoms_step(atoms::ase.atoms, model,walk_params,E_max) #atoms::ase.atoms, model,walk_params,E_max)
    if walk_params.atom_algorithm == "GMC"
        return ase.do_GMC!(atoms,model,walk_params,E_max)
    elseif walk_params.atom_algorithm == "MC"
        return ase.do_MC!(atoms,model,walk_params,E_max)
    elseif walk_params.atom_algorithm == "MD"
        return ase.do_MD!(atoms,model,walk_params,E_cutoff = E_max,KE_cutoff = walk_params.KE_max)
    end
    
end

function do_cell_shear_step(atoms::ase.atoms, model,walk_params,E_max;execute = true,check_energy = true)
   # println("Doing cell shear")
    (p,T) = propose_shear_step(atoms,walk_params.shear_step_size)
    if do_cell_step!(atoms,p,T,E_max,walk_params,model,execute = execute,check_energy=check_energy)
        return (1,1)
    else
        return (1,0)
    end
end

function do_cell_stretch_step(atoms::ase.atoms, model,walk_params,E_max;execute=true,check_energy = true)
   # println("Doing cell stretch")
    (p,T) = propose_stretch_step(atoms,walk_params.stretch_step_size)
    if do_cell_step!(atoms,p,T,E_max,walk_params,model,execute = execute,check_energy=check_energy)
        return (1,1)
    else
        return (1,0)
    end
end


function do_cell_volume_step(atoms::ase.atoms, model,walk_params,E_max;execute=true,check_energy = true)
    #println("Doing cell volume")
    (p,T) = propose_volume_step(atoms,walk_params.volume_step_size)
    if do_cell_step!(atoms,p,T,E_max,walk_params,model,execute = execute,check_energy = check_energy)
        return (1,1)
    else
        return (1,0)
    end
end

function do_cell_step!(atoms,p_accept,T,E_max,ns_params,model;execute = true,check_energy = true)
    nAtoms = length(atoms.positions)
    if (p_accept < 1.0) & (rand() > p_accept)
#        println("rejected for probability")
        return false   
    end

    newCell = atoms.lVecs * T
    newVolume = ase.cell_volume(atoms.latpar * newCell)
#    println("Proposed Volume: ", newVolume)
#    println("Current Volume", ase.cell_volume(atoms))
#    println("Volume limit: ",ns_params.max_volume_per_atom *nAtoms)
    if newVolume > ns_params.max_volume_per_atom * nAtoms
#        println("Volume too large")
#        println("new volume: ", newVolume)
#        println("old volume: ", ase.cell_volume(atoms.latpar * atoms.lVecs))
#        println("Max Volume: ", ns_params.max_volume_per_atom * nAtoms)
        return false
    end
    if ase.min_aspect_ratio(newCell) < ns_params.min_aspect_ratio

 #       println("Aspect ratio too small")
 #       println("new aspect ratio: ", ase.min_aspect_ratio(newCell))
 #       println("old aspect ratio: ", ase.min_aspect_ratio(atoms.lVecs))
 #       println("Max aspect ratio: ", ns_params.min_aspect_ratio)
        return false
    end

    if check_energy
        oldCell = atoms.lVecs
        oldPositions = atoms.positions
   #     println("before")
    #    display(ase.eval_energy(atoms,model))
    
        atoms.lVecs = newCell
        ase.set_cell!(atoms,newCell,scale_atoms = true)
 #       println("after ")
 #       display(ase.eval_energy(atoms,model))
        newEnergy = ase.eval_energy(atoms,model,force_recalc = true)
        if newEnergy < E_max
  #          println("Accepted")
            # For tuning the step sizes we want to check to see if the move will be accepted 
            # but not actually do it.
            if execute
                atoms.model_energy = newEnergy
            else
                ase.set_cell!(atoms,oldCell)
            end
            return true
        else
#            println("Rejected for energy")
#            println("New Energy:", newEnergy)
#            println("E-max: ", E_max)
            ase.set_cell!(atoms,oldCell)
            atoms.positions .= oldPositions # Can I just set cell with scale_atoms = true again?
            return false
        end
    else
        # Not checking for energy, passed the volume and aspect ratio test, so it's a good step
        if execute
            ase.set_cell!(atoms,newCell,scale_atoms = true)
        end
        return true
    end

end


function propose_volume_step(atoms::ase.atoms, stepsize::Float64)
    nAtoms = length(atoms.positions)
    dV = rand(Normal(0,stepsize * nAtoms))
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
    p_accept = minimum([1.0, (new_V/orig_V)^nAtoms])
    
    
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



function do_cell_shape_walk!(atoms::ase.atoms, walker_params)
    possibilities =[do_cell_shear_step,do_cell_stretch_step, do_cell_volume_step]
    # I don't need an energy model to walk the cell shape, but the routine I want to use to do it
    # requires it, so I'm just initializing an empty model.
    empty_model = LennardJones.initialize_empty_model()
    #n_try = 0
    #n_accept = 0
    n_try = walker_params.n_single_walker_steps
    acceptance_rates = Dict("stretch"=>0//n_try,"shear"=>0//n_try,"volume"=>0//n_try)
    for i =1:walker_params.n_single_walker_steps
        shuffle!(possibilities)
        for possible_move in possibilities
            #println(walker_params)
            
            (tried,accepted) = possible_move(atoms,empty_model,walker_params,1000.0,check_energy = false)
            acceptance_rates[string(possible_move)[9:end-5]] += accepted//n_try
#            n_accept += accepted
#            n_try += tried
            #            if (p_accept < 1.0) & (rand() > p_accept)
#                atoms.lVecs = transform * atoms.lVecs
#                ase.DirectToCartesian!(atoms)
           #     println("Accepted")
#            end
        end
    end
    println("tried: ", walker_params.n_single_walker_steps, "Accepted: ",acceptance_rates)
end








end