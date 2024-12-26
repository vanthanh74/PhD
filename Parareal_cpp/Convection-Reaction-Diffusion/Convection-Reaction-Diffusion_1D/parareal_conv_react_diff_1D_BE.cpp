#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;

// -----------------------------------------------------------------
// 1D Convection-Reaction-Diffusion
//      u_t = kappa*(u_xx + u_yy) - 2*p1*u_x - 2*p2*u_y + p3*u + f    in [0,L]x[0,L] x (0,T),
//      u(0,y,t) = u(L,y,t)  = 0                 0 <= y <= L, t >= 0
//      u(x,0,t) = u(x,L,t)  = 0                 0 <= x <= L, t >= 0
//      u(x,y,0) = u0(x,y)                   (x,y) in [0,L]x[0,L]
// -----------------------------------------------------------------

// INPUT PARAMETERS:
//   * L: Length of the domain
//   * Nx: Number of spatial points
//   * T: Total simulation time
//   * Nt_fine: Number of fine time steps
//   * Nt_coarse: Number of coarse time steps
//   * max_iter: Maximum Parareal iterations
//   * alpha: Diffusion coefficient
//   * beta: Convection coefficient
//   * gamma: Reaction coefficient
// OUTPUT:
//   * parareal solution and fine sequential solution
//   * L2-Norm Error between them

// -----------------------------------------------------------------  


// Function prototypes
void backwardEuler1D(const Vec &initial, Vec &solution, double dx, double dt, double alpha, double beta, double gamma, int nt, const Vec &f);
void parareal1D(const Vec &initial, Vec &parareal_solution, double dx, double dt, double alpha, double beta, double gamma, int nt, const Vec &f, int max_iter, Vec &errors);
double computeL2Norm1D(const Vec &v1, const Vec &v2);

int main() {
    // Problem setup
    double L = 1.0; // Length of the domain
    int nx = 100;   // Number of spatial points
    double dx = L / (nx - 1);
    double T = 0.1; // Total time
    int nt = 1000;   // Number of time steps
    double dt = T / nt;

    double alpha = 0.01; // Diffusion coefficient
    double beta = 1.0;   // Convection coefficient
    double gamma = 0.1;  // Reaction coefficient

    Vec x(nx);
    Vec f(nx, 1.0); // Source term, f(x) = 1.0
    Vec initial(nx, 0.0);

    for (int i = 0; i < nx; ++i) {
        x[i] = i * dx;
    }

    // Fine solution using backward Euler
    Vec fine_solution(nx);
    backwardEuler1D(initial, fine_solution, dx, dt, alpha, beta, gamma, nt, f);

    // Parareal solution
    Vec parareal_solution(nx);
    Vec errors;
    int max_iter = 20; // Maximum iterations for Parareal
    parareal1D(initial, parareal_solution, dx, dt, alpha, beta, gamma, nt, f, max_iter, errors);

    // Output results for plotting
    ofstream fine_file("fine_solution_1d.txt"), parareal_file("parareal_solution_1d.txt"), error_file("errors_1d.txt");

    for (int i = 0; i < nx; ++i) {
        fine_file << x[i] << " " << fine_solution[i] << endl;
        parareal_file << x[i] << " " << parareal_solution[i] << endl;
    }

    for (int i = 0; i < errors.size(); ++i) {
        error_file << i + 1 << " " << errors[i] << endl;
    }

    fine_file.close();
    parareal_file.close();
    error_file.close();

    cout << "Results written to fine_solution_1d.txt, parareal_solution_1d.txt, and errors_1d.txt." << endl;

    return 0;
}

void backwardEuler1D(const Vec &initial, Vec &solution, double dx, double dt, double alpha, double beta, double gamma, int nt, const Vec &f) {
    int nx = initial.size();
    solution = initial;

    Vec temp(nx);
    double r_diff = alpha * dt / (dx * dx);
    double r_conv = beta * dt / dx;
    double r_react = gamma * dt;

    for (int t = 0; t < nt; ++t) {
        temp = solution;
        solution[0] = 0.0; // Dirichlet boundary condition
        solution[nx - 1] = 0.0;
        for (int i = 1; i < nx - 1; ++i) {
            solution[i] = (temp[i] + r_diff * (temp[i - 1] - 2.0 * temp[i] + temp[i + 1]) - r_conv * (temp[i] - temp[i - 1]) + dt * f[i]) / (1 + r_react);
        }
    }
}

void parareal1D(const Vec &initial, Vec &parareal_solution, double dx, double dt, double alpha, double beta, double gamma, int nt, const Vec &f, int max_iter, Vec &errors) {
    int nx = initial.size();
    Vec coarse_solution(nx);
    Vec fine_solution(nx);

    Vec current(nx), correction(nx);
    parareal_solution = initial;

    for (int k = 0; k < max_iter; ++k) {
        // Coarse solver
        backwardEuler1D(parareal_solution, coarse_solution, dx, dt * 10, alpha, beta, gamma, 1, f);

        // Fine solver
        backwardEuler1D(parareal_solution, fine_solution, dx, dt, alpha, beta, gamma, nt, f);

        // Correction
        for (int i = 0; i < nx; ++i) {
            correction[i] = fine_solution[i] - coarse_solution[i];
        }

        for (int i = 0; i < nx; ++i) {
            parareal_solution[i] += correction[i];
        }

        double error = computeL2Norm1D(fine_solution, parareal_solution);
        errors.push_back(error);

        if (error < 1e-6) {
            break;
        }
    }
}

double computeL2Norm1D(const Vec &v1, const Vec &v2) {
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return sqrt(sum / v1.size());
}
