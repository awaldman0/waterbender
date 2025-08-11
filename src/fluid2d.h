#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

struct Fluid2D {
    int nx = 0, ny = 0;
    double Lx = 1.0, Ly = 1.0, dx = 0.0, dy = 0.0, dt = 0.016;
    std::vector<std::vector<double>> S0, S1, Ssrc, U0x, U0y, U1x, U1y, Fx, Fy;
    double viscosity = 1e-3, diffusion = 1e-4, dissipation = 1e-2;

    void init(int nx_, int ny_, double Lx_, double Ly_, double dt_);
    void clearSources();
    void injectBox(int cx, int cy, int r, double dens, double fx, double fy);
    void step();
    void sampleVelocityWorld(double xw, double yw, double& vx, double& vy) const;
    void worldToGrid(double xw, double yw, double& ix, double& iy) const;

private:
    static void allocate(std::vector<std::vector<double>>& f, int nx, int ny, double v = 0.0);
    static double bilerp(const std::vector<std::vector<double>>& S, double x, double y, int nx, int ny);
    static void addForce(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src, double dt);
    static void linSolve(std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& x0,
        double a, double c, int nx, int ny, int iters = 40);
    static void diffuse(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src,
        double k, double dt, double dx, int nx, int ny);
    static void dissipate(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src,
        double a, double dt, int nx, int ny);
    static void transport(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src,
        const std::vector<std::vector<double>>& Ux, const std::vector<std::vector<double>>& Uy,
        double dt, double dx, double dy, int nx, int ny);
    void project();
};
