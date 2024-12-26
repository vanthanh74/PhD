import numpy as np
import matplotlib.pyplot as plt

# Load data
fine = np.loadtxt("fine_solution.txt")
parareal = np.loadtxt("parareal_solution.txt")
errors = np.loadtxt("errors.txt")



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


# Plot fine and parareal solutions
fig, ax = plt.subplots(figsize=(9, 6))
plt.plot(fine[:, 0], fine[:, 1], label="Fine Solution", linewidth=2)
plt.plot(parareal[:, 0], parareal[:, 1], label="Parareal Solution", linestyle='--')
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.title("Fine vs Parareal Solution")
plt.savefig("Compare_solution_plot.png")

# Plot convergence curve
fig, ax = plt.subplots(figsize=(9, 6))
plt.semilogy(errors[:, 0], errors[:, 1], marker='o', label="L2 Error")
plt.xlabel("Iteration")
plt.ylabel("L2 Error")
plt.title("Convergence of Parareal Algorithm")
plt.grid()
plt.legend()
plt.savefig("convergence_curve.png")
plt.show()
