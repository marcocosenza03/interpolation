#include "Interpolation/interpolation.hh"
#include <iostream>

int main()
{
    Interpolation::Chebyshev::StandardGrid grid(4);

    for(size_t j=0; j < grid._betaj.size(); j++){
        std::cout << grid._betaj[j] << std::endl;
    }

    std::cout << "\n" << std::endl;

    for(size_t j=0; j < grid._tj.size(); j++){
        std::cout << grid._tj[j] << std::endl;
    }

    return 0;
}