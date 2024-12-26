
import numpy as np
import matplotlib as mpl

from ODE_solver_tools import ODE_solver
from plot_tools import generate_plots
from error_tools import compute_parareal_error, compute_error_f_K

# Number of parareal iterations
K=5
# Refinement params
refine = False
rFactor = [int(2 ** (k+1)) for k in np.arange(K)] # [2, 2^2, ..., 2^K]
# Time params
T = float(12.)
Ng = int(60); dt_g = float(T/Ng) # Ng, dt_g: Number of subintervals, time step subinterval
Nf_K = int(rFactor[-2])
dt_f_K =float(dt_g/Nf_K) # Nf_K = 2^K: Number of steps fine solver in subinerval in last parareal it K-1.
Ne   = 100*Nf_K; dt_e = float(dt_g/Ne)  # Ne, dt_e: Number of steps sequential solver, time step subinterval

Nf = int(0); dt_f = 0.                  # Nf, dt_f: Number of steps of the fine solver in a subinterval, time step fine
if refine:
    print("With refinement")
    Nf = rFactor[0]; dt_f = float(dt_g/Nf)
else:
    print("No refinement")
    Nf = Nf_K; dt_f = dt_f_K

print("K = "+str(K))
print("dt_g = "+str(dt_g))
print("dt_f = "+str(dt_f))
print("dt_e = "+str(dt_e))
print("Ng = " + str(Ng))
print("Nf = " + str(Nf))
print("Nf_K = " + str(Nf_K))
print("Ne = " + str(Ne))

# Initial value
init_state = [0., 1.]
# Define solvers
G = ODE_solver(init_state=init_state,dt=dt_g,num_scheme="RK4")
F = ODE_solver(init_state=init_state,dt=dt_f,num_scheme="RK4")
F_K = ODE_solver(init_state=init_state,dt=dt_f_K,num_scheme="RK4")
Exact = ODE_solver(init_state=init_state,dt=dt_e,num_scheme="RK4")

# Sequential solution (exact)
sol = init_state
seq_ref = [ sol ]
for n in range(Ng):
    sol = Exact.propagate(sol, n*dt_e, Ne)
    seq_ref.append(sol)

# Sequential solution (fine parareal it K)
sol = init_state
seq_f_K = [ sol ]
for n in range(Ng):
    sol = F_K.propagate(sol, n*dt_e, Nf_K)
    seq_f_K.append(sol)

# Parareal solution
X  = []
Fk = []
Gk = []
for k in range(K):
    Xk  = [ np.asarray(init_state) ]
    sol = init_state
    if k == 0:
        # Propagate coarse
        for n in range(Ng):
            sol = G.propagate(sol, n * dt_g, int(1))
            Xk.append(sol)
        X.append(Xk)
        print("k = " + str(k) + ", dt_f = " + str(F.dt) + ", Nf = " + str(Nf))
    else:
        Fk = []
        Gk = []
        # Propagate fine and coarse from X[k-1][n]
        for n in range(Ng):
            # Fine
            sol = F.propagate(X[k-1][n], n*dt_g, Nf)
            Fk.append(sol)
            # Coarse
            sol = G.propagate(X[k-1][n], n*dt_g, int(1))
            Gk.append(sol)
        # Parareal sequence
        for n in range(Ng):
            sol =  G.propagate(Xk[n], n*dt_g, int(1)) + Fk[n]-Gk[n]
            Xk.append(sol)
        X.append(Xk)
        # Refine time step in F
        if refine:
            Nf = int(rFactor[k])
            F.dt = float(dt_g/Nf)
        print("k = "+str(k)+", dt_f = "+str(F.dt)+", Nf = "+str(Nf))

# Error f_K vs reference
err_f_K = compute_error_f_K(seq_ref, seq_f_K)
# Error parareal sequence
err_parareal = compute_parareal_error(seq_ref, X, 1)

# Plots
t = np.arange(0.,T+1e-9,dt_g) # 1e-9 to enforce that T is the final time (problem due to floating point overflow)
generate_plots(t, seq_ref, X, err_parareal, err_f_K, refine)
