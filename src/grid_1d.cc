#include "Interpolation/grid_1d.hh"
#include "Interpolation/gauss_kronrod.hh"
#include <set>


// u=a <=> t=+1
// u=b <=> t=-1
double from_ab_to_m1p1(double u ,double a, double b)
{
    return -2 * (u - a) / (b - a) + 1;
}

double from_m1p1_to_ab(double t, double a, double b)
{
    return 0.5 * (1 - t) * (b - a) + a;
}

double from_ab_to_m1p1_der(double a, double b)
{
    return -2 / (b - a);
}

namespace Interpolation 
{
    SingleDiscretizationInfo::SingleDiscretizationInfo(std::vector<double> inter, std::vector<size_t> g_size,
                            std::function<double(double)> to_i_space,
                            std::function<double(double)> to_i_space_der,
                            std::function<double(double)> to_p_space,
                            std::function<double(double)> to_p_space_der) :
        intervals(inter.size() - 1, {0, 0}), intervals_phys(inter.size() - 1, {0, 0}),
        grid_sizes(g_size), to_inter_space(to_i_space), to_inter_space_der(to_i_space_der),
        to_phys_space(to_p_space), to_phys_space_der(to_p_space_der)
    {
        for(size_t i = 0; i < inter.size() - 1; i++){
            intervals_phys[i] = {inter[i], inter[i + 1]};
            intervals[i] = {to_i_space(inter[i]), to_i_space(inter[i + 1])};
        }
    }
    
    Grid1D::Grid1D(const SingleDiscretizationInfo &d_info) : _d_info(d_info) 
    {   
        size = 0;
        if (d_info.intervals.size() == 0) {
            size_li   = 0;
            c_size    = 0;
            c_size_li = 0;
            return;
        }

        // _stored_grids
        std::vector<size_t> grid_sizes_unique = d_info.grid_sizes;
        std::sort(grid_sizes_unique.begin(), grid_sizes_unique.end());
        grid_sizes_unique.erase(std::unique(grid_sizes_unique.begin(), grid_sizes_unique.end()), grid_sizes_unique.end());
        for(size_t i = 0; i < grid_sizes_unique.size(); i++){
            // Chebyshev::StandardGrid chebyshevpoints(grid_sizes_unique[i]);
            // _stored_grids.insert({grid_sizes_unique[i], chebyshevpoints});
            _stored_grids.insert({grid_sizes_unique[i], Chebyshev::StandardGrid(grid_sizes_unique[i])});
        }


        for(size_t i = 0; i < d_info.intervals.size(); i++) {
            size += d_info.grid_sizes[i] + 1;
        }
        size_li   = static_cast<index_t>(size);
        c_size    = size - (d_info.intervals.size() - 1);
        c_size_li = static_cast<index_t>(c_size);

        _weights.resize(size);
        _weights_der.resize(size);
        _weights_sub.resize(size);
        _from_iw_to_ic.resize(size);
        _from_idx_to_inter.resize(size);

        _coord.resize(c_size);
        _coord_inter.resize(c_size);
        _delim_indexes.resize(d_info.intervals.size() + 1);

        _der_matrix.resize(size);
        for(size_t i = 0; i < _der_matrix.size(); i++) {
            _der_matrix[i].resize(size);
        }


        index_t index       = 0; // Indice dei pesi
        index_t index_coord = 0; // Indice delle coordinate
        for (size_t a = 0; a < d_info.intervals.size(); a++) {

            const Chebyshev::StandardGrid &sg = _stored_grids.at(d_info.grid_sizes[a]);
            _delim_indexes[a] = index;

            for (size_t j = 0; j <= d_info.grid_sizes[a]; j++) {
                // NOTE: _coord stores only the UNIQUES coordinates

                if (j != d_info.grid_sizes[a] || (a == d_info.intervals.size() - 1)) {
                    _coord_inter[index_coord] = from_m1p1_to_ab(sg.t(j), d_info.intervals[a].first, d_info.intervals[a].second);
                    _coord[index_coord] = d_info.to_phys_space(_coord_inter[index_coord]);
                }

                _from_iw_to_ic[index]     = index_coord;
                _from_idx_to_inter[index] = a;
                if (j != d_info.grid_sizes[a]) index_coord++;

                // Note: in order to have proper assignement operator
                // I cannot use [this] in the lambdas. Hence, I need to
                // pass the needed variables by value to the lambdas.
                size_t cached_int_size = _d_info.intervals.size();
                std::pair<double, double> cached_inter = _d_info.intervals[a];

                // ------------------------------
                _weights[index] = [a, j, cached_int_size, cached_inter](double u, const Chebyshev::StandardGrid &sg) -> double {
                    double res = 0;

                    bool condition_x = (a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second);
                    if (u >= cached_inter.first && condition_x) {
                        res += sg.poly_weight(from_ab_to_m1p1(u, cached_inter.first, cached_inter.second), j);
                    }

                    return res;
                };

                // ------------------------------
                _weights_der[index] = [a, j, cached_int_size, cached_inter](double u, const Chebyshev::StandardGrid &sg) -> double {
                    double res = 0;

                    bool condition_x = (a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second);
                    if (u >= cached_inter.first && condition_x) {
                    const double dl_dx = from_ab_to_m1p1_der(cached_inter.first, cached_inter.second);
                    const double dw_dl = sg.poly_weight_der(from_ab_to_m1p1(u, cached_inter.first, cached_inter.second), j);
                    res += dw_dl * dl_dx;
                    }

                    return res;
                };
                _weights_sub[index] = [a, j, cached_int_size, cached_inter](double u, const Chebyshev::StandardGrid &sg) -> double {
                    double res = 0;

                    bool condition_x = (a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second);
                    if (u >= cached_inter.first && condition_x) {
                        res += sg.poly_weight(from_ab_to_m1p1(u, cached_inter.first, cached_inter.second), j) - 1;
                    }

                    return res;
                };

                index_t inner_index = 0;
                for (size_t b = 0; b < d_info.intervals.size(); b++) {
                    for (size_t k = 0; k <= d_info.grid_sizes[b]; k++) {
                        _der_matrix[index][inner_index++] = get_der_matrix(a, j, b, k);
                    }
                }
                
                index++;
            }
        }
        _delim_indexes[d_info.intervals.size()] = index;

        _from_ic_to_iw.resize(c_size);
        for (size_t i = 0; i < _from_iw_to_ic.size(); i++) {
            _from_ic_to_iw[_from_iw_to_ic[i]].push_back(i);
        }








        
        /*size_t N = 1; // Number of unique points
        for(size_t i = 0; i < d_info.grid_sizes.size() - 1; i++){
            N += d_info.grid_sizes[i];
        }
        N += d_info.grid_sizes[d_info.grid_sizes.size() - 1] + 1;

        //_coord
        _coord.resize(N - 1, 0);
        for(size_t i = 0; i < d_info.intervals_phys.size(); i++){}*/


        _integral_weights.resize(c_size_li, 0.);

        size_t i_w;

        std::function<double(double)> full_integrand = [&](double u) -> double {
            return d_info.to_phys_space_der(u) * _weights[i_w](u, get_std_grid(i_w));
        };

        for(size_t j = 0; j < size; j++) {
            size_t j_c = _from_iw_to_ic[j];

            i_w = j;

            auto [vmin, vmax] = get_support_weight_aj(j);

            if(std::abs(vmax - vmin) < 1.0e-15) continue;
            _integral_weights[j_c] += GaussKronrod<GK_21>::integrate(full_integrand, vmin, vmax, 1.0e-10, 1.0e-10);
        }

    }

    double Grid1D::get_der_matrix(size_t a, size_t j, size_t b, size_t k) const
    {
        if(a != b) return 0; 
        return _stored_grids.at(_d_info.grid_sizes[a])._Dij[j][k];
    }


} // namespace Interpolation
