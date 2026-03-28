#include "Interpolation/chebyshev_grid.hh"
//#include <stdexce>

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

    for(size_t i=0; i <= p; i++){
        for(size_t j=0; j <= p; j++){
            if(j==i) continue;
            
            _Dij[i][j] = - (_betaj[i] / _betaj[j]) / (_tj[i] - _tj[j]);
        }
    }
}
/*double StandardGrid::poli_weight(double t, size_t j) const 
{
    double den = 0;
    for(size_t i = 0; i <= _p; i++){
        den = den + _betaj[i] / (t - _tj[i]);
    };
    return _betaj[j] / den;
};
double StandardGrid::poli_weight(double t, size_t j, double den) const 
{
    return _betaj[j] / den;
};*/
double StandardGrid::poli_weight(double t, size_t j, double den) const 
{
    if(std::abs(t - _tj[j]) < 1.0e-15) return 1.; 
    double res = 0.;
    res = _betaj[j] / (t - _tj[j]) / den;
    return res;
}
double StandardGrid::poli_weight(double t, size_t j) const 
{
    if(std::abs(t - _tj[j]) < 1.0e-15) return 1.;
    double den = 0.;
    for(size_t j=0;  j <= _p; j++){
        //if(t == _tj[j]) return fj[j + start];
        if(std::abs(t - _tj[j]) < 1.0e-15) return 0.;
        den += _betaj[j] / (t - _tj[j]);
    }
    double res = 0.;
    res = _betaj[j] / (t - _tj[j]) / den;
    return res;
}
double StandardGrid::poli_weight_der(double t, size_t j) const
{
    double res = 0.;
    for(size_t i = 0; i <= _p; i++){
        if(std::abs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
        res += _Dij[j][i] * poli_weight(t, i);
    }
    return res;
}
double StandardGrid::poli_weight_der(double t, size_t j, double den) const
{
    double res = 0.;
    for(size_t i = 0; i <= _p; i++){
        if(std::abs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
        res += _Dij[j][i] * poli_weight(t, i, den);
    }
    return res;
}
vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
    vector_d fj(_p + 1, 0.);
    for(size_t i = 0; i <= _p; i++){
        fj[i] = fnc(_tj[i]);
    }
    return fj;
}
double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const 
{
    if(t<-1 || t>1){
        throw std::domain_error("StandardGrid::interpolate t must be in [-1,1]");
    }
    if(end - start != _p){
        throw std::domain_error("StandardGrid::interpolate end-start should be = to p");
    }

    double den = 0.;
    for(size_t j=0; j <= _p; j++){
        //if(t == _tj[j]) return fj[j + start];
        if(std::abs(t - _tj[j]) < 1.0e-15) return fj[j + start];
        den += _betaj[j] / (t - _tj[j]);
    }
    double res = 0.;
    for(size_t i=0; i <= _p; i++){
        res += poli_weight(t, i, den) * fj[i + start];
    }
    return res; 
}
double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const
{
    if(t<-1 || t>1){
        throw std::domain_error("StandardGrid::interpolate_der t must be in [-1,1]");
    }
    if(end - start != _p){
        throw std::domain_error("StandardGrid::interpolate_der end-start should be = to p");
    }

    double den = 0.;
    for(size_t j = 0; j <= _p; j++){
        if(std::abs(t - _tj[j]) <= 1.0e-15){
            double res = 0.;
            for(size_t i = 0; i <= _p; i++){
                res += _Dij[i][j] * fj[i + start];
            }
            return res;
        }
        
        den += _betaj[j] / (t - _tj[j]);
    }

    double res = 0.;
    for(size_t i = 0; i <= _p; i++){
        res += poli_weight_der(t, i , den) * fj[i + start];
    }
    return res;
}

//Alternative implementation of interpolate_der(...):
/*double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const
{
    if(t<-1 || t>1){
        throw std::domain_error("StandardGrid::interpolate_der t must be in [-1,1]");
    }
    if(end - start != _p){
        throw std::domain_error("StandardGrid::interpolate_der end-start should be = to p");
    }

    vector_d fjtilde(_p + 1, 0.);
    for(size_t i = 0; i <= _p; i++){
        for(size_t j = 0; j <= _p; j++){
            fjtilde[i] += fj[j + start] * _Dij[j][i];
        }
    }
    
    
    double den = 0.;
    for(size_t j=0; j <= _p; j++){
        //if(t == _tj[j]) return fj[j + start];
        if(std::abs(t - _tj[j]) < 1.0e-15) return fjtilde[j];
        den += _betaj[j] / (t - _tj[j]);
    }
    double res = 0.;
    for(size_t i=0; i <= _p; i++){
        res += poli_weight(t, i, den) * fjtilde[i];
    }
    return res; 
}*/

void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
{
    if(end - start != _p){
        throw std::domain_error("StandardGrid::interpolate end-start should be = to p");
    }

    vector_d temp(end - start + 1, 0.);
    for(size_t j = 0; j < _p; j++){
        for(size_t k = start; k <= end; k++){
            temp[j] += fj[k] * _Dij[k][j];
        }
    }

    for(size_t j = start; j <= end; j++){
        fj[j] = temp[j - start];
    }

}

} // namespace Chebyshev
} // namespace Interpolation