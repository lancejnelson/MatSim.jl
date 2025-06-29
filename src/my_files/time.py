import matplotlib.pyplot as plt

with open(r'C:\Users\rexja\OneDrive - BYU-Idaho\Classes\Nelson Research Group\MatSim.jl\src\my_files\NS.out', 'r') as file:
    lines = file.readlines()

# Process the data
x = []
y = []
for i, line in enumerate(lines):
    if i % 2 == 1 and i > 4:
        parts = line.split(" ")
        if len(parts) == 2:
            x.append(float(parts[0]))
            y.append(float(parts[1]))

# Plot the data
plt.plot(x, y)
plt.xlabel('Epoch')
plt.ylabel('Wall time (s)')
plt.title('NS Wall Time vs Epoch')
plt.show()