#include "fluid2d.h"

static inline double clampd(double v, double lo, double hi) { return std::max(lo, std::min(v, hi)); }
static inline int    clampi(int v, int lo, int hi) { return std::max(lo, std::min(v, hi)); }

void Fluid2D::allocate(std::vector<std::vector<double>>& f, int nx, int ny, double v) {
    f.assign(nx, std::vector<double>(ny, v));
}
double Fluid2D::bilerp(const std::vector<std::vector<double>>& S, double x, double y, int nx, int ny) {
    int i0 = (int)x, j0 = (int)y; int i1 = std::min(i0 + 1, nx - 1), j1 = std::min(j0 + 1, ny - 1);
    double sx = x - i0, sy = y - j0;
    double a = (1.0 - sx) * S[i0][j0] + sx * S[i1][j0];
    double b = (1.0 - sx) * S[i0][j1] + sx * S[i1][j1];
    return (1.0 - sy) * a + sy * b;
}
void Fluid2D::addForce(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src, double dt) {
    int nx = (int)dst.size(), ny = (int)dst[0].size();
    for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) dst[i][j] += dt * src[i][j];
}
void Fluid2D::linSolve(std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& x0,
    double a, double c, int nx, int ny, int iters) {
    for (int k = 0; k < iters; k++) {
        for (int i = 1; i < nx - 1; i++) for (int j = 1; j < ny - 1; j++)
            x[i][j] = (x0[i][j] + a * (x[i + 1][j] + x[i - 1][j] + x[i][j + 1] + x[i][j - 1])) / c;
        for (int i = 0; i < nx; i++) { x[i][0] = x[i][1]; x[i][ny - 1] = x[i][ny - 2]; }
        for (int j = 0; j < ny; j++) { x[0][j] = x[1][j]; x[nx - 1][j] = x[nx - 2][j]; }
    }
}
void Fluid2D::diffuse(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src,
    double k, double dt, double dx, int nx, int ny) {
    dst = src; double a = dt * k / (dx * dx);
    linSolve(dst, src, a, 1.0 + 4.0 * a, nx, ny, 40);
}
void Fluid2D::dissipate(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src,
    double a, double dt, int nx, int ny) {
    for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) dst[i][j] = src[i][j] / (1.0 + dt * a);
}
void Fluid2D::transport(std::vector<std::vector<double>>& dst, const std::vector<std::vector<double>>& src,
    const std::vector<std::vector<double>>& Ux, const std::vector<std::vector<double>>& Uy,
    double dt, double dx, double dy, int nx, int ny) {
    for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) {
        double x = i + 0.5, y = j + 0.5;
        double vx = Ux[i][j], vy = Uy[i][j];
        double midx = x - 0.5 * dt * (vx / dx), midy = y - 0.5 * dt * (vy / dy);
        midx = clampd(midx, 0.0, nx - 1.001); midy = clampd(midy, 0.0, ny - 1.001);
        double vx2 = bilerp(Ux, midx, midy, nx, ny), vy2 = bilerp(Uy, midx, midy, nx, ny);
        double x0 = x - dt * (vx2 / dx), y0 = y - dt * (vy2 / dy);
        x0 = clampd(x0, 0.0, nx - 1.001); y0 = clampd(y0, 0.0, ny - 1.001);
        dst[i][j] = bilerp(src, x0, y0, nx, ny);
    }
}
void Fluid2D::project() {
    std::vector<std::vector<double>> div, p; allocate(div, nx, ny, 0.0); allocate(p, nx, ny, 0.0);
    for (int i = 1; i < nx - 1; i++) for (int j = 1; j < ny - 1; j++) {
        double dxu = (U0x[i + 1][j] - U0x[i - 1][j]) / (2.0 * dx);
        double dyv = (U0y[i][j + 1] - U0y[i][j - 1]) / (2.0 * dy);
        div[i][j] = dxu + dyv;
    }
    linSolve(p, div, 1.0, 4.0, nx, ny, 60);
    for (int i = 1; i < nx - 1; i++) for (int j = 1; j < ny - 1; j++) {
        U1x[i][j] = U0x[i][j] - 0.5 * (p[i + 1][j] - p[i - 1][j]) / dx;
        U1y[i][j] = U0y[i][j] - 0.5 * (p[i][j + 1] - p[i][j - 1]) / dy;
    }
    for (int i = 0; i < nx; i++) { U1x[i][0] = U1x[i][1]; U1x[i][ny - 1] = U1x[i][ny - 2]; U1y[i][0] = 0.0; U1y[i][ny - 1] = 0.0; }
    for (int j = 0; j < ny; j++) { U1y[0][j] = U1y[1][j]; U1y[nx - 1][j] = U1y[nx - 2][j]; U1x[0][j] = 0.0; U1x[nx - 1][j] = 0.0; }
}
void Fluid2D::init(int nx_, int ny_, double Lx_, double Ly_, double dt_) {
    nx = nx_; ny = ny_; Lx = Lx_; Ly = Ly_; dt = dt_; dx = Lx / nx; dy = Ly / ny;
    allocate(S0, nx, ny); allocate(S1, nx, ny); allocate(Ssrc, nx, ny);
    allocate(U0x, nx, ny); allocate(U0y, nx, ny);
    allocate(U1x, nx, ny); allocate(U1y, nx, ny);
    allocate(Fx, nx, ny);  allocate(Fy, nx, ny);
}
void Fluid2D::clearSources() {
    for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) { Ssrc[i][j] = 0.0; Fx[i][j] = 0.0; Fy[i][j] = 0.0; }
}
void Fluid2D::injectBox(int cx, int cy, int r, double dens, double fx, double fy) {
    int i0 = clampi(cx - r, 0, nx - 1), i1 = clampi(cx + r, 0, nx - 1);
    int j0 = clampi(cy - r, 0, ny - 1), j1 = clampi(cy + r, 0, ny - 1);
    for (int i = i0; i <= i1; i++) for (int j = j0; j <= j1; j++) { Ssrc[i][j] += dens; Fx[i][j] += fx; Fy[i][j] += fy; }
}
void Fluid2D::worldToGrid(double xw, double yw, double& ix, double& iy) const {
    double u = (xw + 0.5 * Lx) / Lx, v = (yw + 0.5 * Ly) / Ly; ix = clampd(u * nx, 0.0, nx - 1.001); iy = clampd(v * ny, 0.0, ny - 1.001);
}
void Fluid2D::sampleVelocityWorld(double xw, double yw, double& vx, double& vy) const {
    double ix, iy; worldToGrid(xw, yw, ix, iy); vx = bilerp(U1x, ix, iy, nx, ny); vy = bilerp(U1y, ix, iy, nx, ny);
}
void Fluid2D::step() {
    std::swap(U0x, U1x); std::swap(U0y, U1y); std::swap(S0, S1);
    addForce(U0x, Fx, dt); addForce(U0y, Fy, dt);
    transport(U1x, U0x, U0x, U0y, dt, dx, dy, nx, ny);
    transport(U1y, U0y, U0x, U0y, dt, dx, dy, nx, ny);
    diffuse(U0x, U1x, viscosity, dt, dx, nx, ny);
    diffuse(U0y, U1y, viscosity, dt, dx, nx, ny);
    project();
    addForce(S0, Ssrc, dt);
    transport(S1, S0, U1x, U1y, dt, dx, dy, nx, ny);
    diffuse(S0, S1, diffusion, dt, dx, nx, ny);
    dissipate(S1, S0, dissipation, dt, nx, ny);
    clearSources();
}
