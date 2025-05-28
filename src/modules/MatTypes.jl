module MatTypes

using LinearAlgebra
using StaticArrays


currentDir = abspath(@__DIR__)
libPath = relpath(joinpath(currentDir,"../modules"))

# TODO add types here? Not sure this is necessary at the moment

# include(abspath(libPath,"MatTypes.jl"))

end # module MatTypes
