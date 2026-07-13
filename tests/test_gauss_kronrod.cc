#include <iostream>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include "Interpolation/interpolation.hh"

double exp_function(double x) {
    return exp(x);
}

double lin_function(double x) {
    return x;
}

double quad_function(double x) {
    return pow(x, 2);
}

double gauss_function(double x) {
    return exp(-pow(x, 2));
}

double poly_function(double x) {
    double p = x - 2 * pow(x, 2) + 7 * pow(x, 3) + 4 * pow(x, 15);
    return p;
}

double sine_function(double x) {
    return sin(pow(M_PI, 2) * x);
}

/* int main()
{
    Interpolation::GaussKronrod<Interpolation::GK_21>::Eval int1 = Interpolation::GaussKronrod<Interpolation::GK_21>::gauss_kronrod_simplified(exp_function, 0, M_PI);
    printf("%.10f %e \n", int1.res, int1.err);

    Interpolation::GaussKronrod<Interpolation::GK_21>::Eval int2 = Interpolation::GaussKronrod<Interpolation::GK_21>::gauss_kronrod_simplified(quad_function, -1, 1);
    printf("%.10f %e \n", int2.res, int2.err);

    Interpolation::GaussKronrod<Interpolation::GK_21>::Eval poly_int = Interpolation::GaussKronrod<Interpolation::GK_21>::gauss_kronrod_simplified(poly_function, -1, 1);
    printf("poly int: %.10f %e \n", poly_int.res, poly_int.err);

    Interpolation::GaussKronrod<Interpolation::GK_21>::Eval gauss_int = Interpolation::GaussKronrod<Interpolation::GK_21>::gauss_kronrod_simplified(gauss_function, 0, 1);
    printf("gauss int: %.10f %e \n", gauss_int.res, gauss_int.err);

    std::cout << "\n";

    double res = Interpolation::GaussKronrod<Interpolation::GK_21>::integrate(gauss_function, -1, 1, 1.0e-7, 1.0e-15);
    printf("%.15f\n", res);

    double res_rec = Interpolation::GaussKronrod<Interpolation::GK_21>::integrate_rec(gauss_function, -1, 1, 1.0e-7, 1.0e-15);
    printf("%.15f\n", res_rec);


    double res2 = Interpolation::GaussKronrod<Interpolation::GK_21>::integrate(sine_function, 0, M_PI, 1.0e-2, 1.0e-2);
    printf("%.15f\n", res2);

    double res2_rec = Interpolation::GaussKronrod<Interpolation::GK_21>::integrate_rec(sine_function, 0, M_PI, 1.0e-6, 1.0e-6);
    printf("%.15f\n", res2_rec);

    return 0;
} */





void integate_and_print(const std::function<double(double)> &fnc, const std::string &name, const double exact)
{
   using integrator = Interpolation::GaussKronrod<Interpolation::GK_21>;

   const double res1 = integrator::integrate(fnc, 0, 1, 1.0e-10, 1.0e-10);
   const double res2 = integrator::integrate_rec(fnc, 0, 1, 1.0e-10, 1.0e-10);

   std::cout << std::setprecision(10) << std::scientific;
   std::cout << name << "\t";
   std::cout << "Global adaptive: " << res1 - exact << "\t";
   std::cout << "Local adaptive: " << res2 - exact << "\t";
   std::cout << std::endl;
}

int main()
{
   {
      auto fnc = [](double x) {
         return exp(x);
      };
      integate_and_print(fnc, "exp", exp(1) - 1);
   }
   {
      for (size_t i = 0; i < 50; i += 5) {
         auto fnc = [i](double x) {
            return std::pow(x, i);
         };
         integate_and_print(fnc, "x^" + std::to_string(i), 1. / (i + 1.));
      }
   }

   {
      auto fnc = [](double x) {
         return std::cyl_bessel_j(0, 10 * x);
      };
      integate_and_print(fnc, "J_0(10 x) ", 0.10670113039567368);
   }

   {
      auto fnc = [](double x) {
         return std::cyl_bessel_j(0, 30 * x);
      };
      integate_and_print(fnc, "J_0(30 x) ", 0.029474969627515377);
   }

   double res = Interpolation::GaussKronrod<Interpolation::GK_21>::integrate_rec([](double x) {
      return std::cyl_bessel_j(0, 10 * x);
   }, 0, 1, 1.0e-10, 1.0e-10);

   std::cout << "\n";

   std::printf("%.18f\n", res);
   std::printf("0.10670113039567547\n");

   std::cout << "Hello\n";


   return 0;
}
