module enumeration

using LinearAlgebra:Diagonal,LowerTriangular
using Printf
using StaticArrays
#using Crystal
#using LennardJones

struct parent
    title:: String
    bulk:: Bool
    pLV:: SMatrix{3,3,Float64,9}
    nD:: Int64
    dVecs:: Vector{SVector{3,Float64}}
    k :: Int64
    eps:: Float64

end

# deriverated representation of a crystal
struct deriv
    strN:: Int64
    hnfN::Int64
    hnf_degen:: Int64  
    lab_degen:: Int64
    tot_degen:: Int64
    sizeN:: Int64
    n:: Int64
    pgOps:: Int64
    SNF:: SMatrix{3,3,Int64,9}#Diagonal{Int64,Vector{Int64}}
    HNF:: SMatrix{3,3,Int64,9}#LowerTriangular{Int64, Matrix{Int64}}
    L:: SMatrix{3,3,Int64,9}#Matrix{Int64}
    labeling::String
    arrows::String

end




function read_header(file)
#    f = open(file,"r")
#    lines = readlines(f)
#    close(f)
    lines = String[]
    for line in eachline(file)
        push!(lines,string(line))
        if length(lines) >= 15
            break
        end
    end
    title = string(rstrip(split(lines[1],"\n")[1]))
    if lowercase(rstrip(split(lines[2],"\n")[1])) == "bulk"
        bulk = true
    else
        bulk = false
    end 
    bulkorsurf = rstrip(split(lines[2],"\n")[1])
    pLattice = hcat([parse.(Float64,split(split(i,"#")[1])) for i in lines[3:5]]...)  # Store vectors by columns
    nD = parse(Int,rstrip(split(lines[6])[1]))
    dVecs = [parse.(Float64,split(split(i,"#")[1])) for i in lines[7:7+nD-1]]
    k = parse(Int,rstrip(split(lines[7 + nD],"-")[1]))
    eps = parse(Float64,rstrip(split(lines[9 + nD])[1]))

    return parent(title,bulk,pLattice,nD,dVecs,k,eps)

end

function get_single_line(file_path::String, line_number::Int;header = 15)

    open(file_path,"r") do io
        current_line = 0

        for i in 1:line_number + header - 1
            readline(io)
        end
        line = readline(io)
        true_struct = parse(Int64,split(line)[1])
        if true_struct != line_number
            error("Was looking for line number: $line_number and got instead $true_struct .")
        end
        return split(line)
#        while !eof(io)
#            current_line += 1
#            line = readline(io)
#            if current_line == line_number + header
#                return line
#            end
#        end
#        return "nothing"
    end
end

function read_struct(file,strN)
    #struct_line = 0
    #for (idx,line) in enumerate(eachline(file))
    #    if idx > 15 && startswith(lstrip(line),string(strN))
    #        struct_line = line
    #        break
    #    end
    #end
    #if struct_line == 0
    #    error("Didn't find the specified line")
    #    return
    #end

    struct_line = get_single_line(file,strN)
#    f = open(file,"r")
#    lines = readlines(f)
#    close(f)

#    title = rstrip(split(lines[1],"\n")[1])
#    if lowercase(rstrip(split(lines[2],"\n")[1])) == "bulk"
#        bulk = true
#    else
#        bulk = false
#    end 
#    bulkorsurf = rstrip(split(lines[2],"\n")[1])
#    pLattice = hcat([parse.(Float64,split(split(i,"#")[1])) for i in lines[3:5]]...)  # Store vectors by columns
#    nD = parse(Int,rstrip(split(lines[6])[1]))
#    dVecs = [parse.(Float64,split(split(i,"#")[1])) for i in lines[7:7+nD-1]]
#    k = parse(Int,rstrip(split(lines[7 + nD],"-")[1]))
#    eps = parse(Float64,rstrip(split(lines[9 + nD])[1]))
    
    lineNum = strN + 15
    hnfN = parse(Int,struct_line[2])
    hnf_degen = parse(Int,struct_line[3])
    label_degen = parse(Int,struct_line[4])
    total_degen = parse(Int,struct_line[5])
    sizeN = parse(Int,struct_line[6])
    n = parse(Int,struct_line[7])
    pgOps = parse(Int,struct_line[8])
    @views SNF = SMatrix{3,3,Int64,9}(Diagonal(parse.(Int,struct_line[9:11])))
    a = parse(Int,struct_line[12])
    b = parse(Int,struct_line[13])
    c = parse(Int,struct_line[14])
    d = parse(Int,struct_line[15])
    e = parse(Int,struct_line[16])
    f = parse(Int,struct_line[17])
    HNF = SMatrix{3,3,Int64,9}([a 0 0
                                b c 0 
                                d e f])

    @views l = SVector{9,Int64}(parse(Int,s) for s in struct_line[18:26])
    @views l1 = l[1:3]
    @views l2 = l[4:6]
    @views l3 = l[7:9]
    lTransform = vcat(l1', l2', l3')  # Reshape the vector into a 3x3 matrix
    #    l = parse.(Int,struct_line[18:26])
#    lTransform = hcat((l[i:i+2] for i=1:3:7)...)'
    labeling = struct_line[27]
    arrows = try struct_line[28] catch y repeat("0",length(labeling)) end
    return deriv(strN,hnfN,hnf_degen,label_degen,total_degen,sizeN,n,pgOps,SNF,HNF,lTransform,labeling,arrows)
end

end


