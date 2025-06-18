
using LegendrePolynomials

#Evaluate the spherical harmonics.
function Ylm(l,m,θ,ϕ)
    normalizing = sqrt( (2l + 1)/(4π) * factorial(l - m)/factorial(l + m) )
    return normalizing * Plm(cos(θ),l,m) * exp(1im*m*ϕ)

end

Ylm(1,1,π/2,π/4)



function bispectrum(crystal,cutoff)
    loopBounds = convert.(Int64,cld.(cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] ))
    
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtomType) in enumerate(crystal.atomicBasis), centerAtom in centerAtomType 
            centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
            θ = atan(centerAtomC[2]/centerAtom[1])
            ϕ = atan(centerAtom[3]/ centerAtom)
            ljvals .+= singleAtomLJ(crystal,centerAtomC,iCenter,loopBounds,params,cutoff)    # These two loops are to loop over all possible neighbor atoms
        
    end

end