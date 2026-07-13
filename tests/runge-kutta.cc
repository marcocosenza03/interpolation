#include <iostream>
#include <vector>
#include <functional>


static constexpr double k = 2;

using v2d = std::array<double, 2>;

v2d rhs(double t, v2d q) {
    (void)t;
    return {q[1], -k * q[0]}; // Oscillatore armonico
}

v2d step_rg(double t, double dt, v2d q) {
    v2d k1 = rhs(t, q);

    v2d tmp{};

    tmp[0] = q[0] + k1[0] * dt * 0.5;
    tmp[1] = q[1] + k1[1] * dt * 0.5;
    v2d k2 = rhs(t + dt * 0.5, tmp);

    tmp[0] = q[0] + k2[0] * dt * 0.5;
    tmp[1] = q[1] + k2[1] * dt * 0.5;
    v2d k3 = rhs(t + dt * 0.5, tmp);

    tmp[0] = q[0] + k3[0] * dt;
    tmp[1] = q[1] + k3[1] * dt;
    v2d k4 = rhs(t + dt, tmp);

    /* double k2 = rhs(t + dt * 0.5, q + k1 * dt * 0.5);
    double k3 = rhs(t + dt * 0.5, q + k2 * dt * 0.5);
    double k4 = rhs(t + dt, x + k3 * dt); */
    tmp[0] = q[0] + (dt / 6.) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    tmp[1] = q[1] + (dt / 6.) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

    return tmp;
}

/* double euler(double t, double dt, double x) {
    return x + dt * rhs(t, x);
} */

v2d step_e(double t, double dt, v2d q) {
    v2d k1 = rhs(t, q);
    v2d tmp{};

    tmp[0] = q[0] + k1[0] * dt;
    tmp[1] = q[1] + k1[1] * dt;

    return tmp;
}

int main() 
{
    v2d x0 = {1., 0.};
    double t0 = 0.;
    double tf = 13.;
    size_t nstep = 10000;
    double dt = (tf - t0) / nstep;

    v2d x_rg = x0;
    v2d x_e = x0;
    for(size_t i = 1; i <= nstep; i++) {
        x_rg = step_rg(t0 + i * dt, dt, x_rg);
        x_e = step_e(t0 + i * dt, dt, x_e);
    }

    printf("Runge-Kutta:\n");
    printf("Posizione: %f\nVelocità: %f\n", x_rg[0], x_rg[1]);
    printf("\n");
    printf("Euler:\n");
    printf("Posizione: %f\nVelocità: %f\n", x_e[0], x_e[1]);



    return 0;
}