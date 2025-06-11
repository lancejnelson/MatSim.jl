using LinearAlgebra

A = rand(2000, 2000)
B = rand(2000, 2000)

# Precompile the matrix multiplication
# @time A*B

# Single thread
# begin
#     BLAS.set_num_threads(1)
#     @show BLAS.get_num_threads()
#     @time A*B
# end

# # All threads on the machine
# begin
#     BLAS.set_num_threads(Sys.CPU_THREADS)
#     @show BLAS.get_num_threads()
#     @time A*B
# end


function slow_function(x)
    #BLAS.set_num_threads(Sys.CPU_THREADS)
    
    #@time A*B
    return sqrt(x)^1.5 + log(x + 1)^2 
end

function compute_sum_serial(arr)
    total = 0.0
    for x in arr
        total += slow_function(x)
    end
    return total
end

arr = rand(100000000) .* 100.0  

@time compute_sum_serial(arr)

using Base.Threads

function compute_sum_parallel(arr)
    total = zeros(Float64, nthreads())  # one total per thread
    @threads for i in eachindex(arr)
        total[threadid()] += sqrt(arr[i])^1.5 + log(arr[i] + 1)^2
    end
    return sum(total)
end

@time compute_sum_parallel(arr)
@show nthreads()



