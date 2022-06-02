#ifndef TR7_TR7MATH_H
#define TR7_TR7MATH_H
#pragma once


#include <type_traits>
#include <cmath>
#include <complex>
#include <limits>
#include <algorithm>
#include <numeric>
#include <array>

#include "tr7/tr7defs.h"
#include "tr7/tr7AD.h"


//
// Defines some short & useful math functions and constants.
//


namespace tr7
{
    //
    // Typical numeric constants.
    //

    template<typename T>
    inline constexpr_variable T pi_T = T{ 3.1415926535897932385L };

#ifdef _MSC_VER
    // Visual Studio 2015 has a bug involving variable templates!! It sometimes
    // sets them to zero... A workaround is to perform an explicit instantiation
    template const float pi_T<float>;
    template const double pi_T<double>;
    template const long double pi_T<long double>;

#endif

    // handy names when using doubles, our usual case
    inline constexpr_variable double pi = pi_T<double>;


    //
    // Math functions.
    //

    template<typename T,
             typename TolType = T>
    constexpr_ pure_function bool eq_fp(T x, T y, TolType tol = TolType{ 1 })
    {
        //
        // Kopriva's AlmostEqual routine,
        // tests equality of two floating-point
        // numbers. Input tol is a multiple of
        // "std::numeric_limits<T>::epsilon()".
        //

        const auto teps = tol * std::numeric_limits<TolType>::epsilon();
        const auto diff = AD_abs(x - y);
        const auto ax = AD_abs(x);
        const auto ay = AD_abs(y);

        if (ax <= teps || ay <= teps) {
            return diff <= TolType{ 2 } * teps;
        }
        else {
            return diff <= teps * ax && diff <= teps * ay;
        }

    }


    template<typename T>
    constexpr pure_function T sign(T x)
    {
        //
        // Strict sign function: sign(0) = 0;
        // if a > 0, sign(a) = 1; elseif a < 0, sign(a) = -1.
        //

        return T{ (x > T{ 0 }) - (x < T{ 0 }) };

    }

    template<typename T>
    constexpr pure_function bool samesign(T x, T y)
    {
        //
        // Returns true if x and y have the same sign (0 considered positive).
        //

        return (x >= T{ 0 }) == (y >= T{ 0 });
    }

    template<typename T>
    constexpr pure_function bool samesign(T x, T y, T z)
    {
        //
        // Returns true if x, y and z have the same sign (0 considered positive).
        //

        const auto sx = x >= T{ 0 };

        return (sx == (y >= T{ 0 })) && (sx == (z >= T{ 0 }));
    }


    template<typename T>
    constexpr_ pure_function T mod_mat(T x, T y)
    {
        //
        // Matlab's "mod" function, should only be called for floating-point variables.
        //

        // NEW VERSION FROM MATLAB CODE GENERATION
        auto ret{ x };
        if (std::isfinite(x) && std::isfinite(y)) {
            if (eq_fp(x, T{ 0 })) {
                ret = T{ 0 };
            }
            else {
                if (!eq_fp(y, T{ 0 })) {
                    ret = std::fmod(x, y);
                    bool rEQ0{ eq_fp(ret, T{ 0 }) };
                    if ((!rEQ0) && (y > std::floor(y))) {
                        auto q{ std::abs(x / y) };
                        rEQ0 = (std::abs(q - std::floor(q + T{ 0.5 })) <= T{ 2.2204460492503131E-16 } * q); // 2.220446049250313E-16 == eps, but we'll leave it like this since Matlab hardcodes it
                    }

                    if (rEQ0) {
                        ret = T{ 0 };
                    }
                    else {
                        if ((x < T{ 0 }) != (y < T{ 0 })) {
                            ret += y;
                        }
                    }
                }
            }
        }
        else {
            if (!eq_fp(y, T{ 0 })) {
                ret = std::numeric_limits<T>::quiet_NaN();
            }
        }

        return ret;

    }


    template<typename T>
    constexpr_ pure_function T wrap_to_pi(T ang)
    {
        //
        // Wraps an angle (rad) to [-pi, pi).
        //

        //return mod_mat(ang + pi_T<T>, T{ 2.0 } * pi_T<T>) - pi_T<T>;

#ifdef TR7_WITH_AD
        if constexpr (!is_AD_scalar_v<T>) {
#endif
            return ang - std::floor(ang / (T{ 2 } * pi_T<T>) + T{ 0.5 }) * (T{ 2 } * pi_T<T>);

#ifdef TR7_WITH_AD
        }
        else {
            // with AD, we want the argument of the "floor"
            // function to be of a non-AD-type, while the
            // input "ang" remains an AD-type (the "angle-wrap"
            // operation has derivative = 1 by definition, so
            // doing this won't affect our results)
            return ang - std::floor(AD_force_scalar_value(ang / (T{ 2 } * pi_T<T>) + T{ 0.5 })) * (T{ 2 } * pi_T<T>);

        }
#endif

    }

    template<typename T>
    constexpr_ pure_function T wrap_to_2pi(T ang)
    {
        //
        // Wraps an angle (rad) to [0, 2 * pi).
        //

        return mod_mat(ang, T{ 2 } * pi_T<T>);
    }

    template<typename T>
    constexpr_ pure_function T wrap_to_180(T ang)
    {
        //
        // Wraps an angle (deg) to [-180, 180).
        //

        return mod_mat(ang + T{ 180 }, T{ 360 }) - T{ 180 };
    }

    template<typename T>
    constexpr_ pure_function T wrap_to_360(T ang)
    {
        //
        // Wraps an angle (deg) to [0, 360).
        //

        return mod_mat(ang, T{ 360 });
    }


    template<typename T>
    constexpr pure_function bool intervals_do_not_overlap(T min_a, T max_a, T min_b, T max_b)
    {
        //
        // Returns true if 1-D intervals [min_a, max_a] and [min_b, max_b] don't overlap.
        //

        return (min_a > max_b) || (min_b > max_a);

    }

    template<typename T>
    constexpr pure_function bool intervals_overlap(T min_a, T max_a, T min_b, T max_b)
    {
        //
        // Returns true if 1-D intervals [min_a, max_a] and [min_b, max_b] do overlap.
        //

        return (min_a <= max_b) && (min_b <= max_a);

    }

    template<typename T>
    constexpr pure_function T distance_between_intervals(T min_a, T max_a, T min_b, T max_b)
    {
        //
        // Returns the distance between two 1-D intervals. If positive, they do not overlap.
        //

        return std::max(min_a, min_b) - std::min(max_a, max_b);

    }


    template<typename T,
             typename Arr = std::array<T, 2>>
    constexpr_ std::pair<Arr, bool> sin_cos_solve(T lhs_s, T lhs_c, T rhs)
    {
        //
        // Solves "lhs_s * sin(x) + lhs_c * cos(x) = rhs" for x (2 solutions),
        // only if the equation has real solutions. Returns both solutions in
        // an array and an extra flag that is true if they're valid (AND real),
        // false otherwise (in that case, the array will contain [nan, nan]).
        //

        if (std::abs(lhs_s) > std::numeric_limits<T>::epsilon() && std::abs(lhs_c) > std::numeric_limits<T>::epsilon()) {
            // lhs_s * sin(x) + lhs_c * cos(x) = d * sin(x + p) = rhs
            // d = sqrt(lhs_s^2 + lhs_c^2)
            // p = atan2(lhs_c, lhs_s)
            const auto arg{ rhs / std::sqrt(lhs_s * lhs_s + lhs_c * lhs_c) };
            if (std::abs(arg) <= 1.) {
                const auto p{ std::atan2(lhs_c, lhs_s) };
                const auto x1 = std::asin(arg);

                return std::make_pair(Arr{ wrap_to_pi(x1 - p), wrap_to_pi(pi_T<T> - x1 - p) },
                                      true);

            }

        }
        else if (std::abs(lhs_s) < std::numeric_limits<T>::epsilon() && std::abs(lhs_c) > std::numeric_limits<T>::epsilon()) {
            // lhs_c * cos(x) = rhs
            const auto arg = rhs / lhs_c;
            if (std::abs(arg) <= 1.) {
                const auto x1 = std::acos(arg);

                return std::make_pair(Arr{ wrap_to_pi(x1), wrap_to_pi(-x1) },
                                      true);

            }

        }
        else if (std::abs(lhs_s) > std::numeric_limits<T>::epsilon() && std::abs(lhs_c) < std::numeric_limits<T>::epsilon()) {
            // lhs_s * sin(x) = rhs
            const auto arg = rhs / lhs_s;
            if (std::abs(arg) <= 1.) {
                const auto x1 = std::asin(arg);

                return std::make_pair(Arr{ wrap_to_pi(x1), wrap_to_pi(pi_T<T> - x1) },
                                      true);

            }

        }


        // if we've arrived here, the equation was either "0 = rhs"
        // or had complex roots... invalid solutions
        return std::make_pair(Arr{ std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN() }, false);

    }


    //
    // Conversions between classic units, they
    // don't justify another header...
    //

    template<typename T>
    constexpr pure_function T rad2deg(T ang)
    {
        constexpr_variable auto fac = T{ 180. / pi_T<T> };
        return ang * fac;
    }

    template<typename T>
    constexpr pure_function T deg2rad(T ang_deg)
    {
        constexpr_variable auto fac = T{ pi_T<T> / 180. };
        return ang_deg * fac;
    }


    template<typename T>
    constexpr pure_function T ft2m(T len_ft)
    {
        constexpr_variable auto fac = T{ 0.3048 };
        return len_ft * fac;
    }

    template<typename T>
    constexpr pure_function T m2ft(T len_m)
    {
        constexpr_variable auto fac = T{ 1. / 0.3048 };
        return len_m * fac;
    }


    template<typename T>
    constexpr pure_function T kn2mps(T vel_kn)
    {
        constexpr_variable auto fac = T{ 0.514444 };
        return vel_kn * fac;
    }

    template<typename T>
    constexpr pure_function T mps2kn(T vel_mps)
    {
        constexpr_variable auto fac = T{ 1. / 0.514444 };
        return vel_mps * fac;
    }


    template<typename T>
    constexpr pure_function T kg2N(T mass_kg)
    {
        constexpr_variable auto fac = T{ 9.80665 };
        return mass_kg * fac;
    }

    template<typename T>
    constexpr pure_function T N2kg(T weight_N)
    {
        constexpr_variable auto fac = T{ 1. / 9.80665 };
        return weight_N * fac;
    }


}


#endif
