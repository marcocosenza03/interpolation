#include "Interpolation/interpolation.hh"
#include "Interpolation/ran2.hh"
#include <iostream>


struct DYGrid
{
    DYGrid(size_t nOgata, size_t nXi) : n0(nOgata), nXi(nXi), data(nOgata * nXi, 0.) {}

    double &operator()(size_t i, size_t j)
    {
        return data[i * nXi + j];
    }

    const double &operator()(size_t i, size_t j) const
    {
        return data[i * nXi + j];
    }

    size_t n0, nXi;
    std::vector<double> data;
};

struct Table
{
    std::vector<double> qT; // Experimental points
    std::vector<double> xi_grid;
    std::vector<double> ogata_grid;

    std::vector<DYGrid> WDY;
};

Table compute(const std::vector<double> &qT, std::pair<double, double> bin, size_t nOgata)
{
    using namespace Interpolation;

    if(bin.first >= bin.second) {};

    std::vector<double> bounds;
    std::vector<size_t> sizes;
    
    Grid1D grid();
    Table table;
    table.qT = qT;
    // table.xi_grid = grid._coord;
}


int main() 
{




    return 0;
}