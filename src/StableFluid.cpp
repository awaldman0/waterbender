#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <string>


template <typename T>
T clamp(T val, T minVal, T maxVal) {
    return std::max(minVal, std::min(val, maxVal));
}

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


const int NDIM = 2; // Use 3 for 3D

// Grid parameters
int N[NDIM];       // Number of cells in each dimension
double L[NDIM];    // Physical length of domain
double D[NDIM];    // Cell size = L[i] / N[i]
double O[NDIM];    // Origin (typically 0)
double dt;         // Time step

// Scalar fields (e.g., density)
std::vector<std::vector<double>> S0;
std::vector<std::vector<double>> S1;
std::vector<std::vector<double>> Ssource;

std::vector<std::vector<double>> U0[NDIM]; // velocity field
std::vector<std::vector<double>> U1[NDIM]; // velocity field'
std::vector<std::vector<double>> F[NDIM];

// Fluid parameters
double viscosity = 0.001;
double diffusionconstant = 0.0001;
double dissipationrate = 0.01;




// Allocation helper for 2D fields
void allocateField(std::vector<std::vector<double>>& field, int nx, int ny) {
    field.resize(nx, std::vector<double>(ny, 0.0));
}

// Clamp helper
template <typename T>
T clamp(T val, T minVal, T maxVal) {
    return std::max(minVal, std::min(val, maxVal));
}

// Bilinear interpolation
double bilerp(const std::vector<std::vector<double>>& S, double x, double y) {
    int i0 = (int)x;
    int j0 = (int)y;
    int i1 = std::min(i0 + 1, N[0] - 1);
    int j1 = std::min(j0 + 1, N[1] - 1);
    double sx = x - i0;
    double sy = y - j0;

    double a = (1 - sx) * S[i0][j0] + sx * S[i1][j0];
    double b = (1 - sx) * S[i0][j1] + sx * S[i1][j1];
    return (1 - sy) * a + sy * b;
}

// RK2 backtrace
void TraceParticle(double x, double y,
    std::vector<std::vector<double>> U[NDIM],
    double dt,
    double& x0, double& y0) {
    double vx1 = U[0][(int)x][(int)y];
    double vy1 = U[1][(int)x][(int)y];

    double midx = x - 0.5 * dt * vx1 / D[0];
    double midy = y - 0.5 * dt * vy1 / D[1];

    midx = clamp(midx, 0.0, N[0] - 1.001);
    midy = clamp(midy, 0.0, N[1] - 1.001);

    double vx2 = U[0][(int)midx][(int)midy];
    double vy2 = U[1][(int)midx][(int)midy];

    x0 = x - dt * vx2 / D[0];
    y0 = y - dt * vy2 / D[1];

    x0 = clamp(x0, 0.0, N[0] - 1.001);
    y0 = clamp(y0, 0.0, N[1] - 1.001);
}

// Transport
void Transport(std::vector<std::vector<double>>& dst,
    const std::vector<std::vector<double>>& src,
    std::vector<std::vector<double>> U[NDIM],
    double dt) {
    for (int i = 0; i < N[0]; ++i) {
        for (int j = 0; j < N[1]; ++j) {
            double x = i + 0.5;
            double y = j + 0.5;
            double x0, y0;
            TraceParticle(x, y, U, dt, x0, y0);
            dst[i][j] = bilerp(src, x0, y0);
        }
    }
}

// Add force
void addForce(std::vector<std::vector<double>>& field,
    const std::vector<std::vector<double>>& force,
    double dt) {
    for (int i = 0; i < N[0]; ++i)
        for (int j = 0; j < N[1]; ++j)
            field[i][j] += dt * force[i][j];
}

// Linear solver
void linSolve(std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& x0,
    double a, double c, int iterations = 20) {
    for (int k = 0; k < iterations; ++k) {
        for (int i = 1; i < N[0] - 1; ++i) {
            for (int j = 1; j < N[1] - 1; ++j) {
                x[i][j] = (x0[i][j] +
                    a * (x[i + 1][j] + x[i - 1][j] +
                        x[i][j + 1] + x[i][j - 1])) / c;
            }
        }
    }
}

// Diffuse
void Diffuse(std::vector<std::vector<double>>& dst,
    const std::vector<std::vector<double>>& src,
    double k, double dt) {
    dst = src; // initial guess
    double a = dt * k / (D[0] * D[0]); // assumes dx = dy
    linSolve(dst, src, a, 1 + 4 * a);
}

// Dissipate
void Dissipate(std::vector<std::vector<double>>& dst,
    const std::vector<std::vector<double>>& src,
    double a, double dt) {
    for (int i = 0; i < N[0]; ++i)
        for (int j = 0; j < N[1]; ++j)
            dst[i][j] = src[i][j] / (1.0 + dt * a);
}

// Project
void Project(std::vector<std::vector<double>> U1[NDIM],
    std::vector<std::vector<double>> U0[NDIM],
    double dt) {
    std::vector<std::vector<double>> div(N[0], std::vector<double>(N[1], 0.0));
    std::vector<std::vector<double>> p(N[0], std::vector<double>(N[1], 0.0));

    // Compute divergence
    for (int i = 1; i < N[0] - 1; ++i) {
        for (int j = 1; j < N[1] - 1; ++j) {
            double dx = (U0[0][i + 1][j] - U0[0][i - 1][j]) / (2 * D[0]);
            double dy = (U0[1][i][j + 1] - U0[1][i][j - 1]) / (2 * D[1]);
            div[i][j] = dx + dy;
        }
    }

    linSolve(p, div, 1.0, 4.0);

    // Subtract gradient
    for (int i = 1; i < N[0] - 1; ++i) {
        for (int j = 1; j < N[1] - 1; ++j) {
            U1[0][i][j] = U0[0][i][j] - 0.5 * (p[i + 1][j] - p[i - 1][j]) / D[0];
            U1[1][i][j] = U0[1][i][j] - 0.5 * (p[i][j + 1] - p[i][j - 1]) / D[1];
        }
    }
}

// Velocity step
void Vstep(std::vector<std::vector<double>> U1[NDIM],
    std::vector<std::vector<double>> U0[NDIM],
    double visc, std::vector<std::vector<double>> F[NDIM],
    double dt) {
    for (int i = 0; i < NDIM; ++i)
        addForce(U0[i], F[i], dt);

    for (int i = 0; i < NDIM; ++i)
        Transport(U1[i], U0[i], U0, dt);

    for (int i = 0; i < NDIM; ++i)
        Diffuse(U0[i], U1[i], visc, dt);

    Project(U1, U0, dt);
}

// Scalar step
void Sstep(std::vector<std::vector<double>>& S1,
    std::vector<std::vector<double>>& S0,
    double k, double a,
    std::vector<std::vector<double>> U[NDIM],
    std::vector<std::vector<double>>& source,
    double dt) {
    addForce(S0, source, dt);
    Transport(S1, S0, U, dt);
    Diffuse(S0, S1, k, dt);
    Dissipate(S1, S0, a, dt);
}

void saveDensityImage(const std::vector<std::vector<double>>& field, const std::string& filename) {
    int width = N[0];
    int height = N[1];
    std::vector<unsigned char> image(width * height);

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            double val = clamp(field[i][j], 0.0, 1.0); // assumes field in [0,1]
            image[(height - 1 - j) * width + i] = static_cast<unsigned char>(val * 255.0);
        }
    }

    stbi_write_png(filename.c_str(), width, height, 1, image.data(), width);
}


int main() {
    // Grid resolution and spacing
    N[0] = 128;
    N[1] = 128;
    L[0] = 1.0;
    L[1] = 1.0;
    D[0] = L[0] / N[0];
    D[1] = L[1] / N[1];
    O[0] = 0.0;
    O[1] = 0.0;
    dt = 0.1;



    // Allocate scalar fields
    allocateField(S0, N[0], N[1]);
    allocateField(S1, N[0], N[1]);
    allocateField(Ssource, N[0], N[1]);

    // Allocate velocity and force fields
    for (int i = 0; i < NDIM; ++i) {
        allocateField(U0[i], N[0], N[1]);
        allocateField(U1[i], N[0], N[1]);
        allocateField(F[i], N[0], N[1]);
    }

    // Main simulation loop
    bool simulating = true;
    int frame = 0;

    while (simulating) {
        // === User-defined forces and sources ===
        // Example: inject density and upward force at center
        int cx = N[0] / 2;
        int cy = N[1] / 2;
        for (int i = -2; i <= 2; ++i) {
            for (int j = -2; j <= 2; ++j) {
                int x = cx + i;
                int y = cy + j;
                if (x >= 0 && x < N[0] && y >= 0 && y < N[1]) {
                    Ssource[x][y] = 5.0;
                    F[1][x][y] = 50.0; // upward force
                }
            }
        }

        // === Swap buffers ===
        for (int i = 0; i < NDIM; ++i)
            std::swap(U1[i], U0[i]);
        std::swap(S1, S0);

        // === Simulation steps ===
        Vstep(U1, U0, viscosity, F, dt);
        Sstep(S1, S0, diffusion, dissipation, U1, Ssource, dt);

        // === Output image every 10 frames ===
        if (frame % 10 == 0) {
            std::string filename = "frame_" + std::to_string(frame) + ".png";
            saveDensityImage(S1, filename);
            std::cout << "Saved " << filename << "\n";
        }

        // === Exit after 100 frames ===
        frame++;
        if (frame > 100) simulating = false;
    }

    return 0;
}


    return 0;
}
