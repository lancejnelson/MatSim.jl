module StillingerWeber

#include("/Users/legoses/OneDrive - BYU-Idaho/codes/Crystal.jl")
#include("/Users/jerem/Code/Computation Material Science/MatSim.jl/src/libraries/ase.jl")
using .CrystalMethods


function isEquivalentAtom(atomOne,atomTwo)  # Atoms are assumed to be in direct coordinates.

    #=If the atoms are separated by a lattice vector we can be sure that they are equivalent atoms (translationally).
       To check this, subtract the two vectors and see if the resulting vector is only composed of 
       integers. 
    =#
    return all(round.(Int64,atomOne - atomTwo) - (atomOne - atomTwo) == [0,0,0]) 


end

function isValidTriplet(triplet,radius)

    r1 = norm(centerAtom - atomTwo)
    r2 = norm(centerAtom - atomThree)

    differentPos =  !all(triplet.centerAtom .≈ triplet.atomTwo) && !all(triplet.centerAtom .≈ triplet.atomThree) && !all(atomTwo .≈ triplet.atomThree) 
    insideRadius = r1 < radius && r2 < radius
    return differentPos && insideRadius

end

function initializeTriplet(centerAtom,atomTwo,atomThree,types)
    diffVecOne = triplet.atomTwo - triplet.centerAtom
    diffVecTwo = triplet.atomThree - triplet.centerAtom
    r1 = norm(diffVecOne)
    r2 = norm(diffVecTwo)
    cosθ = diffVecOne' * diffVecTwo/(norm(diffVecOne) * norm(diffVecTwo))
    return triplet(centerAtom,atomTwo,atomThree,r1,r2,cosθ,types)
    #return [r1,r2,cosθ]

end

#function getRandTheta(triplet)
#    diffVecOne = triplet.atomTwo - triplet.centerAtom
#    diffVecTwo = triplet.atomThree - triplet.centerAtom
#    r1 = norm(diffVecOne)
#    r2 = norm(diffVecTwo)
#    cosθ = diffVecOne' * diffVecTwo/(norm(diffVecOne) * norm(diffVecTwo))
#    return [r1,r2,cosθ]
#
#end


function SW_pair(pair,SW)
    return (SW.A[pair.types...] * pair.r^SW.p[pair.types...] - SW.B[pair.types...] * pair.r^SW.q[pair.types...] ) * exp(SW.δ[pair.types...]/(pair.r - SW.a[pair.types...]))
end

function SW_threeBody(triplet,SW)
    return SW.λ[triplet.types...]* exp(SW.γ[triplet.types[1:2]...]/(triplet.r1 - SW.b[triplet.types[1:2]...])  +  SW.γ[[triplet.types[1],triplet.types[3]]...]/(triplet.r2 - SW.b[[triplet.types[1],triplet.types[3]]...])) * (triplet.cosθ + 1/3)^2
end


# Calculates sum of all pair energies in a crystal.
function pairEnergy(crystal::crystal, params, atomOne, typeTwo, radius)
    convertToDirect!(crystal)
    atom = convertToCartesian(crystal.lVecs,crystal.atomicBasis[atomOne[1]][atomOne[2]])
    pairE = 0.0
    # The three inner loops are to find all of the atoms inside of the cutoff sphere.  The outer loop is
    count = 0
    for atomA in crystal.atomicBasis[typeTwo], i = -3:3, j = -3:3, k = -10:10
        newAtomDirect = atomA + I * [i,j,0]
        newAtomCart = convertToCartesian(crystal.lVecs, newAtomDirect)
        r =  norm(atom - newAtomCart)
        if 0 < r <= radius
            count += 1
            pairE += SW_pair(r,params[string(atomOne[1]) * string(typeTwo)])
        end
    end
    return [count,pairE]

end

function singleAtomEnergy(crystal::crystal,atom,cutoffRadius,model)

    pairE = 0.0
    threeBodyE = 0.0
    pairE += pairEnergy(crystal,model.pairParams,atom,1,cutoffRadius)[2]  # Compute all of the pair interactions between this atom and all other atoms that are of type 1
    pairE += pairEnergy(crystal,model.pairParams,atom,2,cutoffRadius)[2]  # Compute all of the pair interactions between this atom and all other atoms that are of type 2
    threeBodyE += threeBodyEnergy(crystal,model.threeBodyParams,atom,[1,1],cutoffRadius)[2] # Compute all of the three-body interactions between this atom and all other pairs of atoms that are of type 1 and 1
    threeBodyE += threeBodyEnergy(crystal,model.threeBodyParams,atom,[1,2],cutoffRadius)[2] # Compute all of the three-body interactions between this atom and all other pairs of atoms that are of type 1 and 2
    threeBodyE += threeBodyEnergy(crystal,model.threeBodyParams,atom,[2,2],cutoffRadius)[2] # Compute all of the three-body interactions between this atom and all other pairs of atoms that are of type 2 and 2
    return pairE/2.0 + threeBodyE
end

function totalEnergy(crystal::crystal,cutoffRadius,model)
    
    totalEnergy = 0
    for (i,atomType) in enumerate(crystal.atomicBasis)
        for (j,atom) in enumerate(atomType)
            println("Calculating total energy associated with atom:")
            println([i,j])
            addEnergy = singleAtomEnergy(crystal,[i,j],cutoffRadius,model)
            println(addEnergy)
            totalEnergy += addEnergy
        end
    end
    return totalEnergy

end


function threeBodyEnergy(crystal::crystal, params,atomOne, types, radius,allTouches = false)
    convertToDirect!(crystal)  # Put all of the basis vectors in direct coordinates
    atomOneC = convertToCartesian(crystal.lVecs,crystal.atomicBasis[atomOne[1]][atomOne[2]])  # Convert the atom of interest into cartesian coordinates
    atomTwoC = Float64[0,0,0]
    atomThreeC = Float64[0,0,0]
    #r = 0
    #atoms = zeros(3,3)
    #atoms[:,1] .= convertToCartesian(crystal,atom)

    
    perms = [SVector{3,Float64}(n,m,0) for n in -3:3 for m in -3:3] # All possible translates in direct coords.
    if types[1]==types[2]
        allAtoms = [(n,) .+ perms for n in crystal.atomicBasis[types[1]]]  # Get all of the a atom coords (direct)
        allperms = reduce(vcat,allAtoms)  # Combine from all of the various basis vectors into one list
        dubs = combinations(allperms,2)  # Enumerate all possible combinations of two atoms.
    else
        aAtoms = [(n,) .+ perms for n in crystal.atomicBasis[types[1]]] # Get all of the a atom coords (direct)
        bAtoms = [(n,) .+ perms for n in crystal.atomicBasis[types[2]]] # Get all of the b atom coords (direct)
        allA = reduce(vcat,aAtoms) # Combine from all of the various basis vectors into one list
        allB = reduce(vcat,bAtoms) # Combine from all of the various basis vectors into one list
        dubs = Iterators.product(allA,allB)
#        dubs = [[n,m] for n in allA for m in allB] # Enumerate all possible combinations of A-B pairs
    end
    count = 0  # How many triplets did I find
    degenerate = 0  # How many times did I include a triplet that will be double counted in a total energy calculation
    correction = 0.0  # Keep track of the correction factor so that I can subtract if off if I'm doing a total energy calculation.
    pot = 0.0    #  Keep track of the total function value.
      
    
    for pair in dubs
        thisTriplet = getRandTheta(atomOneC,atomTwoC,atomThreeC)
        atomTwoC .= convertToCartesian(crystal.lVecs, pair[1])
        atomThreeC .= convertToCartesian(crystal.lVecs, pair[2])

        if isValidTriplet(atomOneC, atomTwoC,atomThreeC,radius)
            count += 1
            pot += SW_threeBody(getRandTheta(atomOneC,atomTwoC,atomThreeC)...,params[string(atomOne[1])*string(types[1])* string(types[2])])
            # Include this one
        else # If one of the arms of the triplet is too far away, go to the next possibility
            continue
        end 

        # The code below is in case you want to find **every three body** that touches this particular atom... not just the three bodies that have this atom at it's center. If you do this and repeat for all of the other atoms in the unit cell, you would triple count all of the interactions and have to correct for that to get the total energy.  

        # I can envision a scenario where you want to know how the total energy changes when moving a single atom.  Instead of calculating the energy of the whole cell twice (once before the move and once after)  it will be more efficient to just calculate the affect of a single atom moving.

        # If one of the atoms in the three-body is equivalent to the center atom, then we don't need to translate the center of the three-body to that site; You're going to get the exact same interaction eventually.  Otherwise, you need to re-center the three body on the other two atoms and include them.
        if allTouches && !isEquivalentAtom(crystal.atomicBasis[atomOne[1]][atomOne[2]],pair[1])
            degenerate += 1
            siteEnergy = SW_threeBody(getRandTheta(atomTwoC,atomOneC,atomThreeC)...,params[string(types[1])*string(atomOne[1])* string(types[2])])
            correction +=  siteEnergy
            count += 1
            pot += siteEnergy
        end

        if allTouches && !isEquivalentAtom(crystal.atomicBasis[atomOne[1]][atomOne[2]],pair[2])
            degenerate += 1
            siteEnergy = SW_threeBody(getRandTheta(atomThreeC,atomOneC,atomTwoC)...,params[string(types[2])*string(atomOne[1])* string(types[1])])
            correction +=  siteEnergy
            count += 1
            pot += siteEnergy
        end
    end
    return [count,pot,correction,degenerate]
end

end