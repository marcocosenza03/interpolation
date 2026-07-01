#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "Interpolation/zeros.hh"

static constexpr double h = 10e-3;

double psi(double x) 
{
    return x * std::tanh(M_PI * 0.5 * std::sinh(x));
}

double psi_prime(double x)
{
    double num = M_PI * x * std::cosh(x) + std::sinh(M_PI * std::sinh(x));
    double den = 1 + std::cosh(M_PI * std::sinh(x));

    return num / den;
}

std::vector<double> tabulated_nodes()
{
    std::vector<double> nodes(1000, 0.);
    for(size_t i = 0; i < nodes.size(); i++) {
        nodes[i] = (M_PI / h) * psi(h * J0_zero[i] / M_PI);
    }
    return nodes;
}

/* std::vector<double> tabulated_weights()
{
    std::vector<double> weights(1000, 0.);
    for(size_t i = 0; i < weights.size(); i++) {
        double tmp = M_PI * psi_prime(h * J0_zero[i] / M_PI);
        tmp *= std::cyl_bessel_j()
    }
} */

double ogata(const std::function<double(double)> &fnc, const std::vector<double> &nodes ,const std::vector<double> &weights, double eps = 1.0e-5)
{
    double res = 0.;
    for(size_t i = 0; i < nodes.size(); i++) {
        double tmp = fnc(nodes[i]) * weights[i];
        if(i > 10 && std::abs(tmp) / std::abs(res) < eps) break;

        res += tmp;
    }
    return res;
}

int main() 
{
    auto nodes = tabulated_nodes();
    // auto weights = tabulated_weights(nodes);

    /* {
        auto fnc = [](double b) {
            return b * exp(- b * b / 2.) / 2.;
        };
        std::cout << ogata(fnc, nodes, weights, 1.0e-6) << std::endl;
    } */

    // std::cout << std::cyl_bessel_j(0, 2) << "\n";


    return 0;
}