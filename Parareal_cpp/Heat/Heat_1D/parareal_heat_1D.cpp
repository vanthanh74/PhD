#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;


// ---------------------------------------------------------------
// 1D Heat equation
// u_t = alpha*u_xx 			 in [0,L]x[0,T],
// u(x,0) = u0(x)                  in (0,L),
// u(0,t) = u(L,t) = 0            t in (0,T)
// -----------------------------------------------------------------

// INPUT PARAMETERS:
//   * L: Length of the domain
//   * Nx: Number of spatial points
//   * T: Total simulation time
//   * Nt_fine: Number of fine time steps
//   * Nt_coarse: Number of coarse time steps
//   * max_iter: Maximum Parareal iterations
//   * alpha: Diffusion coefficient
// OUTPUT:
//   * parareal solution and fine sequential solution
//   * L2-Norm Error between them
   
// -----------------------------------------------------------------


// Parameters
const double L = 1.0;       // Length of the domain
const int Nx = 100;         // Number of spatial points
const double T = 0.1;       // Total simulation time
const int Nt_fine = 1000;   // Number of fine time steps
const int Nt_coarse = 100;  // Number of coarse time steps
const int max_iter = 20;    // Maximum Parareal iterations
const double alpha = 0.01;  // Diffusion coefficient

// Functions
void initialize(Vec &u) {
    for (int i = 0; i < u.size(); ++i) {
        double x = i * L / (u.size() - 1);
        u[i] = sin(M_PI * x); // Initial condition
    }
}

// Dirichlet boundary condition
void apply_boundary_conditions(Vec &u) {
    u[0] = u[u.size() - 1] = 0.0;
}


// Explicit Forward Euler fine solver
void fine_solver(const Vec &u0, Vec &u, double dt, int steps) {
    int N = u0.size();
    u = u0;
    Vec u_next(N);

    double dx = L / (N - 1);
    double r = alpha * dt / (dx * dx);

    for (int step = 0; step < steps; ++step) {
        for (int i = 1; i < N - 1; ++i) {
            u_next[i] = u[i] + r * (u[i - 1] - 2 * u[i] + u[i + 1]);
        }
        apply_boundary_conditions(u_next);
        u.swap(u_next);
    }
}

// Explicit Forward Euler coarse one-step solver 
void coarse_solver(const Vec &u0, Vec &u, double dt, int steps) {
    fine_solver(u0, u, dt, 1); // Coarse solver is less accurate
}

// Compute L2-norm error 
void compute_error(const Vec &u1, const Vec &u2, double &error) {
    error = 0.0;
    for (int i = 0; i < u1.size(); ++i) {
        error += pow(u1[i] - u2[i], 2);
    }
    error = sqrt(error / u1.size());
}


// Main function
int main() {
    int N = Nx;
    int Nt = Nt_coarse;
    double dx = L / (N - 1);
    double dt_fine = T / Nt_fine;
    double dt_coarse = T / Nt_coarse;

    // Initial condition
    Vec u0(N);
    initialize(u0);
    
    // Fine solver reference solution
    Vec u_fine(N);
    fine_solver(u0, u_fine, dt_fine, Nt_fine);

    // Parareal solution
    Mat u_parareal(Nt + 1, Vec(N));
    u_parareal[0] = u0;

    for (int t = 1; t <= Nt; ++t) {
        coarse_solver(u_parareal[t - 1], u_parareal[t], dt_coarse, 1);
    }

    vector<double> errors(max_iter);
    for (int iter = 0; iter < max_iter; ++iter) {
        for (int t = 1; t <= Nt; ++t) {
            Vec u_fine_temp(N);
            fine_solver(u_parareal[t - 1], u_fine_temp, dt_fine, Nt_fine / Nt);

            Vec u_coarse_temp(N);
            coarse_solver(u_parareal[t - 1], u_coarse_temp, dt_coarse, 1);

            for (int i = 0; i < N; ++i) {
                u_parareal[t][i] += u_fine_temp[i] - u_coarse_temp[i];
            }
        }

        // Compute L2 error
        compute_error(u_fine, u_parareal[Nt], errors[iter]);
    }

    // Write results for plotting
    ofstream fine_file("fine_solution.dat");
    ofstream parareal_file("parareal_solution.dat");
    ofstream error_file("errors.dat");

    for (int i = 0; i < N; ++i) {
        fine_file << i * dx << " " << u_fine[i] << endl;
    }

    for (int i = 0; i < N; ++i) {
        parareal_file << i * dx << " " << u_parareal[Nt][i] << endl;
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        error_file << iter << " " << errors[iter] << endl;
    }

    fine_file.close();
    parareal_file.close();
    error_file.close();

    cout << "Simulation complete. Results written to files." << endl;

    return 0;
}
