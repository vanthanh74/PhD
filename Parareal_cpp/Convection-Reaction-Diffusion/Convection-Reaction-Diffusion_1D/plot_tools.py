import numpy as np
import matplotlib.pyplot as plt

# Load data
fine_solution = np.loadtxt("fine_solution_1d.txt")
parareal_solution = np.loadtxt("parareal_solution_1d.txt")
errors = np.loadtxt("errors_1d.txt")

# Extract data
x_fine = fine_solution[:, 0]
y_fine = fine_solution[:, 1]
x_para = parareal_solution[:, 0]
y_para = parareal_solution[:, 1]

# Plot fine solution
plt.figure()
plt.plot(x_fine, y_fine, label="Fine Solution", linestyle='-', marker='o')
plt.title("Fine Solution")
plt.xlabel("x")
plt.ylabel("Solution")
plt.legend()
plt.grid()
plt.savefig("fine_solution_plot.png")
#plt.show()

# Plot Parareal solution
plt.figure()
plt.plot(x_para, y_para, label="Parareal Solution", linestyle='-', marker='x')
plt.title("Parareal Solution")
plt.xlabel("x")
plt.ylabel("Solution")
plt.legend()
plt.grid()
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
