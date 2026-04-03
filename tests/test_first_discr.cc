#include "Interpolation/interpolation.hh"

double testfnc(double x)
{
    return exp(x);
}

int main()
{
    Interpolation::Chebyshev::StandardGrid grid(3);
    auto fj = grid.discretize(testfnc);
    double approx = grid.interpolate(0.27, fj, 0, fj.size() - 1);
    double exact = testfnc(0.27);

    std::cout << "Zeroth derivative" << std::endl;
    std::cout << approx << std::endl;
    std::cout << exact << std::endl;

    double approxder = grid.interpolate_der(0.27, fj, 0, fj.size() - 1);
    double exactder = testfnc(0.27);

    std::cout << "First derivative" << std::endl;
    std::cout << approxder << std::endl;
    std::cout << exactder << std::endl;

    std::cout << "\n";

    std::cout << "Punti tj:" << std::endl;
    for(size_t i=0; i<grid._tj.size(); i++){
        std::cout << grid._tj[i] << std::endl;
    }

    std::cout << "\n";

    std::cout << "Componenti di fj:" << std::endl;
    for(size_t i=0; i<fj.size(); i++){
        std::cout << fj[i] << std::endl;    
    }

    std::cout << "\n";

    grid.apply_D(fj, 0, fj.size() - 1);

    std::cout << "Componenti di fj dopo che ha agito apply_D(...):" << std::endl;
    for(size_t i=0; i<fj.size(); i++){
        std::cout << fj[i] << std::endl;    
    }

    /*
    std::cout << "Componenti di Dij:" << std::endl;
    for(size_t i = 0; i < fj.size(); i++){
        std::cout << grid._Dij[i][i] << std::endl;
    }
    */
    
    return 0;
}