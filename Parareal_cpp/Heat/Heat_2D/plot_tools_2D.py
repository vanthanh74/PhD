import numpy as np
import matplotlib.pyplot as plt

# Load data
fine_solution = np.loadtxt("fine_solution_2d.txt")
parareal_solution = np.loadtxt("parareal_solution_2d.txt")
errors = np.loadtxt("errors_2d.txt")

# Reshape data for 2D plotting
nx = int(np.sqrt(len(fine_solution)))  # Assuming square grid
ny = nx
x = fine_solution[:, 0].reshape((nx, ny))
y = fine_solution[:, 1].reshape((nx, ny))
fine_values = fine_solution[:, 2].reshape((nx, ny))
parareal_values = parareal_solution[:, 2].reshape((nx, ny))



# Plotting
#axis and title fontsize
params = {'axes.labelsize': 12,
          'axes.titlesize': 12,
          'xtick.labelsize': 9,   # Size of the x tick labels
          'ytick.labelsize': 9    # Size of the y tick labels
}
          
plt.rcParams.update(params)

#legend fontsize
plt.rc('legend',fontsize=12) # using a size in points

# Plot fine solution
fig, ax = plt.subplots(figsize=(9, 6))
plt.contourf(x, y, fine_values, levels=50, cmap='viridis')
plt.colorbar(label="Fine Solution")
plt.title("Fine Solution")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("fine_solution_plot.png")
#plt.show()

# Plot Parareal solution
fig, ax = plt.subplots(figsize=(9, 6))
plt.contourf(x, y, parareal_values, levels=50, cmap='viridis')
plt.colorbar(label="Parareal Solution")
plt.title("Parareal Solution")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("parareal_solution_plot.png")
#plt.show()

# Plot convergence curve
iterations = errors[:, 0]
error_values = errors[:, 1]
fig, ax = plt.subplots(figsize=(9, 6))
plt.semilogy(iterations, error_values, marker='o')
plt.title("Convergence Curve")
plt.xlabel("Iteration")
plt.ylabel("L2 Norm Error")
plt.grid(True)
plt.savefig("convergence_curve.png")
plt.show()
