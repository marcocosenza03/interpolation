#include "Interpolation/interpolation.hh"

// Equispaced points in the interval [-1,1]
double fnc1(size_t n, size_t p)
{
    return -1 + (2 / (p + 1)) * n;
}

double fnc2(size_t n, size_t p)
{
    return cos((p - n) * M_PI / p);
}

double testfunction(double x)
{
    return exp(x);
}


int main()
{
    /*
    Interpolation::Chebyshev::StandardGrid grid(4);

    for(size_t j=0; j < grid._betaj.size(); j++){
        std::cout << grid._betaj[j] << std::endl;
    }

    std::cout << "\n" << std::endl;
    

    for(size_t j=0; j < grid._tj.size(); j++){
        std::cout << grid._tj[j] << std::endl;
    }
    */
    
    
    
    // size_t deg = 7;

    // Interpolation::Chebyshev::StandardGrid grid1(deg);
    /*std::cout << "Chebyshev:" << std::endl;
    for(size_t i = 0; i < grid1._tj.size(); i++){
        std::cout << "_tj["<< i <<"] = " << grid1._tj[i] << std::endl;
    }
    for(size_t i = 0; i < grid1._tj.size(); i++){
        std::cout << "_betaj["<< i <<"] = " << grid1._betaj[i] << std::endl;
    
    for(size_t i = 0; i < grid1._tj.size(); i++){
        for(size_t j = 0; j < grid1._tj.size(); j++){
            std::cout << "_Dij["<< i <<"]["<< j <<"] = " << grid1._Dij[i][j] << std::endl;
        }
    }*/
    
    /*Interpolation::Generic::StandardGrid grid2(fnc2, deg);
    std::cout << "Generic:" << std::endl;
    for(size_t i = 0; i < grid2._tj.size(); i++){
        std::cout << "_tj["<< i <<"] = " << grid2._tj[i] << std::endl;
    }*/
    /*for(size_t i = 0; i < grid2._tj.size(); i++){
        std::cout << "_lambdaj["<< i <<"] = " << grid2._lambdaj[i] << std::endl;
    }
    for(size_t i = 0; i < grid2._tj.size(); i++){
        for(size_t j = 0; j < grid2._tj.size(); j++){
            std::cout << "_Dij["<< i <<"]["<< j <<"] = " << grid2._Dij[i][j] << std::endl;
        }
    }
    for(size_t i = 0; i < grid2._tj.size(); i++){
        std::cout << "_lambdaj["<< i <<"] = " << grid2._lambdaj[i] << std::endl;
    }*/

    /*for(size_t i = 0, j = 1; i <= 10; i++, j++){
        std::cout << " "<< i <<" "<< j <<" " << std::endl;  
    }

    Interpolation::vector_d testfunctionvalues = grid2.discretize(testfunction);
    for(size_t i = 0; i < grid2._tj.size(); i++){
        std::cout << testfunctionvalues[i] << std::endl;
    }*/

    // printf("f(%d)\n", 1); 

    //double value = grid2.interpolate(0.5, testfunctionvalues, 0, grid2._p, Interpolation::Generic::StandardGrid::STRATEGY::SBF);

    //std::cout << value << std::endl; 

    //std::cout << "Front: " << grid2._tj.front() << std::endl; 

    /*Interpolation::vector_d vec;
    vec.resize(deg + 1, 0.);
    for(size_t i = 0; i <= deg; i++){
        vec[i] = cos(i * M_PI / deg);
    }
    
    Interpolation::Generic::StandardGrid grid2(vec);
    std::cout << "Generic:" << std::endl;
    for(size_t i = 0; i < vec.size(); i++){
        std::cout << grid2._tj[i] << std::endl;
    }*/ 

    std::vector<double> points = {1, 2, 3, 4};
    std::vector<size_t> degrees = {4, 2, 5, 5};

    Interpolation::SingleDiscretizationInfo info(points, degrees, testfunction, testfunction, testfunction, testfunction);
    for(size_t i = 0; i < info.grid_sizes.size(); i++){
        std::cout << info.grid_sizes[i] << " ";
    }
    std::cout << "\n";
    for(size_t i = 0; i < info.intervals_phys.size(); i++){
        std::cout << info.intervals_phys[i].first << " " << info.intervals_phys[i].second << std::endl;
    }

    Interpolation::Grid1D grid1d(info);
    std::cout << grid1d._stored_grids.size() << std::endl;
    std::cout << grid1d._stored_grids[4]._p << std::endl;
    // std::cout << grid1d._stored_grids[1]._p << std::endl; 
    /*for(size_t i = 0; i < grid1d._coord.size(); i++){
        std::cout << grid1d._coord[i] << " ";
    }
    std::cout << "\n";*/

    return 0;
}