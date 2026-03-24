#include "Interpolation/chebyshev_grid.hh"

namespace Interpolation
{
namespace Chebyshev
{
StandardGrid::StandardGrid(size_t p)
{
	_p = p;

    for(size_t j=0; j<=p; j++){
        _tj.push_back(cos(j * M_PI / static_cast<double>(p)));
    }

    for(size_t j=0; j<=p; j++){
        double sign = j % 2 == 0 ? +1 : -1;
        double scaling;
        if(j==0 || j==p) scaling = 0.5;
        else scaling = 1.;
        _betaj.push_back(sign * scaling);
    }

    _Dij.resize(p+1, vector_d(p + 1, 0.));
    // D is a vector of dimensions p+1 of vectors of dimensions p+1 of zeros
    _Dij[0][0] = (2. * p * p + 1) / 6.;
    _Dij[p][p] = -_Dij[0][0];

    for(size_t j = 1; j<p; j++){
        _Dij[j][j] = -0.5 * _tj[j] / (1. - pow(_tj[j], 2));
    }

    for(size_t i=0; i<= p; i++){
        for(size_t j=0; i<= p; i++){
            if(j==i) continue;

            _Dij[i][j] = - (_betaj[i] / _betaj[j]) / (_tj[i] - _tj[j]);
        }
    }
}
} // namespace Chebyshev
} // namespace Interpolation