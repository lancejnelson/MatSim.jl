# Read data from NS.out, skipping the first 4 rows, and assign columns to V and E

filename = "NS.out"

# Read the data, skipping the first 4 rows (header)
data = readdlm(filename, skipstart=4)

# Assign columns to V and E, then reverse both arrays to maintain data pairs
V = data[:, 1] # Fractional volume in configuration space
E = data[:, 2] # Energy in eV

kB = 1.0              # Boltzmann constant (set to 1.0 for reduced units)
# kB = 1.380649e-23     # J/K
# kB = 8.617333262145e-5  # eV/K
T = 300.0              # Temperature in K

# Calculate ΔV (difference between adjacent V values)
ΔV = V[1:end-1] - V[2:end]

boltzmann_factors = exp.(-(E / (kB * T)))

# Weighted sum for the partition function
partition_function = sum(ΔV .* boltzmann_factors[1:end-1])
# TODO figure out if it should be boltzmann_factors[2:end] or boltzmann_factors[1:end-1]

using Plots

# plot temperature dependence
T_values = 0:10:3000  # Temperatures from 100 K to 1000 K in steps of 10 K
partition_functions = Float64[]
enthalpies = Float64[]
heat_capacities = Float64[]
for T in T_values
    print("Temperature: ", T, " K\n")
    boltzmann_factors = exp.(-(E / (kB * T)))
    
    partition_function = sum(ΔV .* boltzmann_factors[1:end-1])
    println("Partition function: ", partition_function)
    push!(partition_functions, partition_function)

    # energy or enthalpy from partition function
    enthalpy = -kB * T * log(partition_function)
    println("Enthalpy: ", enthalpy)
    push!(enthalpies, enthalpy)

    # heat capacity from partition function
    heat_capacity = kB * T^2 * (log(partition_function) + 1)
    println("Heat capacity: ", heat_capacity)
    push!(heat_capacities, heat_capacity)
    print("\n")

end
# Plot the results
plot(T_values, partition_functions, xlabel="Temperature (K)", ylabel="Partition Function", title="Partition Function vs Temperature", legend=:topright)
plot(T_values, enthalpies, xlabel="Temperature (K)", ylabel="Enthalpy (eV)?", title="Enthalpy vs Temperature", legend=:topright)
plot(T_values, heat_capacities, xlabel="Temperature (K)", ylabel="Heat Capacity (eV/K)?", title="Heat Capacity vs Temperature", legend=:topright)

# plot enthalpy vs heat capacity
plot(enthalpies, heat_capacities, xlabel="Enthalpy (eV)", ylabel="Heat Capacity (eV/K)", title="Enthalpy vs Heat Capacity", legend=:topright)