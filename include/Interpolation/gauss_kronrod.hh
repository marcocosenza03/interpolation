#pragma once

#include "Interpolation/default.hh"
#include <queue>

namespace Interpolation
{

namespace cpt
{

template <typename Rule, std::size_t N>
concept IsGKRule = requires(Rule r) {
    { r.size } -> std::convertible_to<size_t>;
    requires std::same_as<decltype(r.x_gk), const std::array<double, N>>;
    requires std::same_as<decltype(r.w_g), const std::array<double, (N - 1) / 2>>;
    requires std::same_as<decltype(r.w_gk), const std::array<double, N>>;
};

} // namespace cpt

struct GK_21 {
    constexpr static size_t size = 11;
    const static std::array<double, 11> x_gk;
    const static std::array<double, 5> w_g;
    const static std::array<double, 11> w_gk;
};

struct GK_41 {
    constexpr static size_t size = 21;
    const static std::array<double, 21> x_gk;
    const static std::array<double, 10> w_g;
    const static std::array<double, 21> w_gk;
};

struct GK_61 {
   constexpr static size_t size = 31;
   const static std::array<double, 31> x_gk;
   const static std::array<double, 15> w_g;
   const static std::array<double, 31> w_gk;
};

constexpr double _LOCAL_DBL_MIN_     = 2.2250738585072014e-308;
constexpr double _1e8_LOCAL_DBL_MIN_ = 2.2250738585072014e-300;
constexpr double _LOCAL_DBL_EPSILON_ = 2.2204460492503131e-16;
constexpr double _LOCAL_DBL_MAX_     = 1.7976931348623157e+308;

static inline bool subinterval_too_small(double a1, double a2, double b2)
{
   const double tmp = (1 + _LOCAL_DBL_EPSILON_) * (std::fabs(a2) + _LOCAL_DBL_MIN_);
   return std::fabs(a1) <= tmp && std::fabs(b2) <= tmp;
}

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
struct GaussKronrod {

    using rule_t = Rule;

    /**
    * @brief Integrate the input function over an interval
    *
    * @param fnc The function to be integrated
    * @param a The lower limit of the integral
    * @param b The upper limit of the integral
    * @param tol_rel The maximum relative error accepted
    * @param tol_abs The maximum absolute error accepted
    * @return double
    */
    static double integrate(const std::function<double(double)> &fnc, double a, double b, double tol_rel = 1.0e-6, double tol_abs = 1.0e-10);

    static double integrate_rec(const std::function<double(double)> &fnc, double a, double b, double tol_rel = 1.0e-6, double tol_abs = 1.0e-10); // recursive integration
        
    public: // make public for testing
    struct Eval {
        double res, err;
    };
    struct Item {
        Eval e;
        double low, high;
    };
    struct ByErrAbsMax {
        bool operator()(const Item &lhs, const Item &rhs) const
        {
            return lhs.e.err < rhs.e.err;
        }
    };
    using PQ = std::priority_queue<Item, std::vector<Item>, ByErrAbsMax>;

    static Eval gauss_kronrod_simplified(const std::function<double(double)> &fnc, double a, double b);
};

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
Interpolation::GaussKronrod<Rule>::Eval GaussKronrod<Rule>::gauss_kronrod_simplified(const std::function<double(double)> &fnc, double a, double b) 
{
    double sum_gk = 0;
    double sum_g = 0;
    double M = (b + a) / 2;
    double D = (b - a) / 2;

    for(size_t i = 0; i < Rule::w_gk.size() - 1; i++) {
        sum_gk += D * Rule::w_gk[i] * (fnc(M + D * Rule::x_gk[i]) + fnc(M - D * Rule::x_gk[i]));
    }
    sum_gk += D * Rule::w_gk[Rule::x_gk.size() - 1] * fnc(M + D * Rule::x_gk[Rule::x_gk.size() - 1]);

    for(size_t i = 0; i < Rule::w_g.size(); i++) {
        sum_g += D * Rule::w_g[i] * (fnc(M + D * Rule::x_gk[2 * i + 1]) + fnc(M - D * Rule::x_gk[2 * i + 1]));
    }

    /* Eval result;
    result.res = sum_gk;
    result.err = std::abs(sum_gk - sum_g); */
    Eval result = {.res = sum_gk, .err = std::abs(sum_gk - sum_g)};
    return result;
}

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
double Interpolation::GaussKronrod<Rule>::integrate(const std::function<double(double)> &fnc, double a, double b, double tol_rel, double tol_abs)
{
    Eval seed = gauss_kronrod_simplified(fnc, a, b);
    
    if(seed.err < tol_abs) return seed.res;
    Item seed_item = {.e = seed, .low = a, .high = b};

    PQ pq;
    pq.emplace(seed_item);

    double result = 0.;
    double error = 0.;

    result += seed.res;
    error += seed.err;

    while(!pq.empty()) {
        Item current = pq.top();
        pq.pop();

        double c = (current.high + current.low ) / 2;

        if(subinterval_too_small(current.low, c, current.high)) continue;

        Eval e_low = gauss_kronrod_simplified(fnc, current.low, c);
        Eval e_high = gauss_kronrod_simplified(fnc, c, current.high);

        Item left = {.e = e_low, .low = current.low, .high = c};
        Item right = {.e = e_high, .low = c, .high = current.high};

        pq.emplace(left);
        pq.emplace(right);

        result += e_high.res + e_low.res - current.e.res;
        error += e_high.err + e_low.err - current.e.err;

        double tol = std::max(tol_abs, tol_rel * std::abs(result));
        if(error < tol) break;
    }
    return result;
}

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
double Interpolation::GaussKronrod<Rule>::integrate_rec(const std::function<double(double)> &fnc, double a, double b, double tol_rel, double tol_abs)
{
    Eval seed = gauss_kronrod_simplified(fnc, a , b);
    if(seed.err < tol_abs) return seed.res;
    
    std::vector<Item> stack;
    Item seed_item = {.e = seed, .low = a, .high = b};
    stack.push_back(seed_item);

    double result = 0;
    result += seed.res;

    while(stack.size() != 0)
    {
        Item current = stack.back();
        stack.pop_back();

        // result -= current.e.high;

        double c = (current.high - current.low) / 2;
        if(subinterval_too_small(current.low, c, current.high)) continue;
        else result -= current.e.res;

        Eval e_low = gauss_kronrod_simplified(fnc, current.low, c);
        Eval e_high = gauss_kronrod_simplified(fnc, c, current.high);

        Item left = {.e = e_low, .low = current.low, .high = c};
        Item right = {.e = e_high, .low = c, .high = current.high};

        double tol_low = std::max(tol_abs, tol_rel * std::abs(e_low.res));
        double tol_high = std::max(tol_abs, tol_rel * std::abs(e_high.res));

        if(e_low.err < tol_low) {
            result += e_low.res;
        }
        else {
            result += e_low.res;
            stack.push_back(left);
        }
        if(e_high.err < tol_high) {
            result += e_high.res;
        }
        else {
            result += e_high.res;
            stack.push_back(right);
        }
    }

    return result;
}

} // namespace Interpolation
