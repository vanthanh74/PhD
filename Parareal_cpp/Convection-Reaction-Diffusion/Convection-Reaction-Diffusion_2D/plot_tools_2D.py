import numpy as np
import matplotlib.pyplot as plt

# Load data
fine_solution = np.loadtxt("fine_solution_2d.txt")
parareal_solution = np.loadtxt("parareal_solution_2d.txt")
errors = np.loadtxt("errors_2d.txt")

# Extract coordinates and reshape data for plotting
nx = int(np.sqrt(len(fine_solution)))  # Assuming square grid
ny = nx
x = fine_solution[:, 0].reshape((nx, ny))
y = fine_solution[:, 1].reshape((nx, ny))
fine_values = fine_solution[:, 2].reshape((nx, ny))
parareal_values = parareal_solution[:, 2].reshape((nx, ny))

# Plot fine solution
plt.figure()
plt.contourf(x, y, fine_values, levels=50, cmap='viridis')
plt.colorbar(label="Fine Solution")
plt.title("Fine Solution")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("fine_solution_plot.png")
#plt.show()

# Plot Parareal solution
plt.figure()
plt.contourf(x, y, parareal_values, levels=50, cmap='viridis')
plt.colorbar(label="Parareal Solution")
plt.title("Parareal Solution")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("parareal_solution_plot.png")
#plt.show()

# Plot Convergence Curve
iterations = errors[:, 0]
error_values = errors[:, 1]
plt.figure()
plt.semilogy(iterations, error_values, marker='o', label="L2 Norm Error")
plt.title("Convergence Curve")
plt.xlabel("Iteration")
plt.ylabel("Error (L2 Norm)")
plt.grid(True)
plt.legend()
plt.savefig("convergence_curve.png")
plt.show()
