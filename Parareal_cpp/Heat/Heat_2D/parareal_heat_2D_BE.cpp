#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;

// -----------------------------------------------------------------
// 2D Heat equation
//      u_t = alpha*(u_xx + u_yy) + f    in [0,L]x[0,L] x (0,T),
//      u(0,y,t) = u(L,y,t)  = 0                 0 <= y <= L, t >= 0
//      u(x,0,t) = u(x,L,t)  = 0                 0 <= x <= L, t >= 0
//      u(x,y,0) = u0(x,y)                   (x,y) in [0,L]x[0,L]
// -----------------------------------------------------------------

// INPUT PARAMETERS:
//   * Lx, Ly: number of coarse intervals in x,y
//   * nx, ny: Number of spatial points in x,y
//   * T: Total time
//   * nt: Number of time steps
//   * max_iter: Maximum Parareal iterations
//   * alpha: Diffusion coefficient
// OUTPUT:
//   * parareal solution and fine sequential solution
//   * L2-Norm Error between them

// -----------------------------------------------------------------

// Function prototypes
void backwardEuler2D(const Mat &initial, Mat &solution, double dx, double dy, double dt, double alpha, int nt, const Mat &f);
void parareal2D(const Mat &initial, Mat &parareal_solution, double dx, double dy, double dt, double alpha, int nt, const Mat &f, int max_iter, Vec &errors);
double computeL2Norm2D(const Mat &m1, const Mat &m2);

int main() {
    // Problem setup
    double Lx = 1.0; // Length of the domain in x
    double Ly = 1.0; // Length of the domain in y
    int nx = 100;     // Number of spatial points in x
    int ny = 100;     // Number of spatial points in y
    double dx = Lx / (nx - 1);
    double dy = Ly / (ny - 1);
    double T = 0.1;  // Total time
    int nt = 1000;    // Number of time steps
    double dt = T / nt;
    double alpha = 0.01; // Diffusion coefficient

    Mat x(nx, Vec(ny)), y(nx, Vec(ny));
    Mat f(nx, Vec(ny, 1.0)); // Source term, f(x, y) = 1.0
    Mat initial(nx, Vec(ny, 0.0));

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            x[i][j] = i * dx;
            y[i][j] = j * dy;
        }
    }

    // Fine solution using backward Euler
    Mat fine_solution(nx, Vec(ny));
    backwardEuler2D(initial, fine_solution, dx, dy, dt, alpha, nt, f);

    // Parareal solution
    Mat parareal_solution(nx, Vec(ny));
    Vec errors;
    int max_iter = 20; // Maximum iterations for Parareal
    parareal2D(initial, parareal_solution, dx, dy, dt, alpha, nt, f, max_iter, errors);

    // Output results for plotting
    ofstream fine_file("fine_solution_2d.txt"), parareal_file("parareal_solution_2d.txt"), error_file("errors_2d.txt");

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            fine_file << x[i][j] << " " << y[i][j] << " " << fine_solution[i][j] << endl;
            parareal_file << x[i][j] << " " << y[i][j] << " " << parareal_solution[i][j] << endl;
        }
    }

    for (int i = 0; i < errors.size(); ++i) {
        error_file << i + 1 << " " << errors[i] << endl;
    }

    fine_file.close();
    parareal_file.close();
    error_file.close();

    cout << "Results written to fine_solution_2d.txt, parareal_solution_2d.txt, and errors_2d.txt." << endl;

    return 0;
}

void backwardEuler2D(const Mat &initial, Mat &solution, double dx, double dy, double dt, double alpha, int nt, const Mat &f) {
    int nx = initial.size();
    int ny = initial[0].size();
    solution = initial;

    Mat temp(nx, Vec(ny));
    double rx = alpha * dt / (dx * dx);
    double ry = alpha * dt / (dy * dy);

    for (int t = 0; t < nt; ++t) {
        temp = solution;
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                solution[i][j] = (temp[i][j] + rx * (temp[i - 1][j] - 2.0 * temp[i][j] + temp[i + 1][j]) + ry * (temp[i][j - 1] - 2.0 * temp[i][j] + temp[i][j + 1]) + dt * f[i][j]) / (1 + 2 * (rx + ry));
            }
        }
        // Dirichlet boundary conditions
        for (int i = 0; i < nx; ++i) {
            solution[i][0] = 0.0;
            solution[i][ny - 1] = 0.0;
        }
        for (int j = 0; j < ny; ++j) {
            solution[0][j] = 0.0;
            solution[nx - 1][j] = 0.0;
        }
    }
}

void parareal2D(const Mat &initial, Mat &parareal_solution, double dx, double dy, double dt, double alpha, int nt, const Mat &f, int max_iter, Vec &errors) {
    int nx = initial.size();
    int ny = initial[0].size();
    Mat coarse_solution(nx, Vec(ny));
    Mat fine_solution(nx, Vec(ny));

    Mat current(nx, Vec(ny)), correction(nx, Vec(ny));
    parareal_solution = initial;

    for (int k = 0; k < max_iter; ++k) {
        // Coarse solver
        backwardEuler2D(parareal_solution, coarse_solution, dx, dy, dt * 10, alpha, 1, f);

        // Fine solver
        backwardEuler2D(parareal_solution, fine_solution, dx, dy, dt, alpha, nt, f);

        // Correction
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                correction[i][j] = fine_solution[i][j] - coarse_solution[i][j];
            }
        }

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                parareal_solution[i][j] += correction[i][j];
            }
        }

        double error = computeL2Norm2D(fine_solution, parareal_solution);
        errors.push_back(error);

        if (error < 1e-6) {
            break;
        }
    }
}

double computeL2Norm2D(const Mat &m1, const Mat &m2) {
    double sum = 0.0;
    int nx = m1.size();
    int ny = m1[0].size();
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            sum += (m1[i][j] - m2[i][j]) * (m1[i][j] - m2[i][j]);
        }
    }
    return sqrt(sum / (nx * ny));
}
