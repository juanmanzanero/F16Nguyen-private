#ifndef F16_NGUYEN_LINEARIZE_PLANT_H
#define F16_NGUYEN_LINEARIZE_PLANT_H
#pragma once


#include <vector>

#include "cppad/cppad.hpp"

#include "F16_Nguyen/contiguous_ranges.h"


//
// Defines the "linearize_plant" function,
// that linearizes a plant using the
// "CppAD" library.
//

namespace F16_Nguyen
{
    template<typename T>
    struct linearization_results
    {
        //
        // Holds the results that function
        // "linearize_plant" calculates. All
        // matrices are stored in column-major
        // format.
        //

        using value_type = T;


        std::vector<value_type> A; // state matrix, of size "num_states x num_states"
        std::vector<value_type> B; // input matrix, of size "num_states x num_inputs"
        std::vector<value_type> C; // output matrix, of size "num_outputs x num_states"
        std::vector<value_type> D; // feedthrough matrix, of size "num_outputs x num_inputs"

        std::vector<value_type> xdot; // vector of time derivatives of the states at the linearization point
        std::vector<value_type> y; // vector of outputs at the linearization point
        std::vector<value_type> x; // vector of states that define the linearization point
        std::vector<value_type> u; // vector of inputs that define linearization point

    };


    template<typename ValueTypeOfResults = double,
             class Plant, typename StatesType, typename InputsType>
    inline linearization_results<ValueTypeOfResults> linearize_plant(const Plant &plant,
                                                                     const StatesType &x0, const InputsType &u0)
    {
        //
        // Calculates the state-space representation
        // of a plant's linearization:
        //
        //  "y0   = plant.outputs(x0, u0),
        //  xdot0 = plant.derivatives(x0, y0, u0),
        //
        //  Dx    = x - x0,
        //  Du    = u - u0,
        //  Dy    = y - y0,
        //  Dxdot = xdot - xdot0,
        //
        //  Dxdot = A * Dx + B * Du,
        //  Dy    = C * Dx + D * Du",
        //
        // at a given point "[x0, u0]" using the
        // "CppAD" library, and returns a
        // "linearization_results" object with
        // the results.
        //

        using value_type = ValueTypeOfResults;
        using AD_value_type = CppAD::AD<value_type>;


        // get the sizes of the plant
        const auto num_states = plant.num_states();
        const auto num_outputs = plant.num_outputs();
        const auto num_inputs = plant.num_inputs();


        // ensure sizes
        ensure_minimum_size(num_states, x0);
        ensure_minimum_size(num_inputs, u0);


        // copy the states and inputs into a single vector of
        // AD "independent variables", and declare them as such
        std::vector<AD_value_type> xu0_AD(num_states + num_inputs);
        const auto x0_AD = xu0_AD.begin();
        const auto u0_AD = std::next(xu0_AD.begin(), num_states);

        std::copy_n(x0.begin(), num_states, x0_AD);
        std::copy_n(u0.begin(), num_inputs, u0_AD);

        CppAD::Independent(xu0_AD);


        // evaluate the "dependent variables"  calling the
        // plant's state-space member functions, and put them
        // into a single vector
        std::vector<AD_value_type> yxdot0_AD(num_outputs + num_states);
        const auto y0_AD = yxdot0_AD.begin();
        const auto xdot0_AD = std::next(yxdot0_AD.begin(), num_outputs);

        std::copy_n(plant.outputs(x0_AD, u0_AD).cbegin(), num_outputs, y0_AD);
        std::copy_n(plant.derivatives(x0_AD, y0_AD, u0_AD).cbegin(), num_states, xdot0_AD);


        // create the AD functions and stop the recording
        CppAD::ADFun<value_type> funs_AD;
        funs_AD.Dependent(xu0_AD, yxdot0_AD);


        // we're ready to evaluate the AD functions, now with numerical arguments
        struct CppAD_Value_const_iterator : std::vector<AD_value_type>::const_iterator
        {
            using base_type = typename std::vector<AD_value_type>::const_iterator;
            CppAD_Value_const_iterator(base_type it) : base_type{ it } {}
            auto operator*() const { return CppAD::Value(base_type::operator*()); }
        };

        std::vector<value_type> xu0(CppAD_Value_const_iterator(xu0_AD.cbegin()),
                                    CppAD_Value_const_iterator(xu0_AD.cend()));

        const auto yxdot0 = funs_AD.Forward(0, xu0);
        const auto Jac0 = funs_AD.Jacobian(xu0);


        // separate the resulting row-major jacobian into
        // the "A", "B", "C", and "D" col-major
        // state-space matrices, they're stored like this
        //
        //   "Jac0 = [C, D;
        //            A, B]".
        //
        std::vector<value_type> C(num_outputs * num_states);
        std::size_t count = 0u;
        for (auto i = 0u; i < num_outputs; ++i) {
            for (auto j = 0u; j < num_states; ++j) {
                C[i + j * num_outputs] = Jac0[count++];
            }
            count += num_inputs;

        }

        std::vector<value_type> D(num_outputs * num_inputs);
        count = num_states;
        for (auto i = 0u; i < num_outputs; ++i) {
            for (auto j = 0u; j < num_inputs; ++j) {
                D[i + j * num_outputs] = Jac0[count++];
            }
            count += num_states;

        }

        std::vector<value_type> A(num_states * num_states);
        count = num_outputs * (num_states + num_inputs);
        for (auto i = 0u; i < num_states; ++i) {
            for (auto j = 0u; j < num_states; ++j) {
                A[i + j * num_states] = Jac0[count++];
            }
            count += num_inputs;

        }

        std::vector<value_type> B(num_states * num_inputs);
        count = num_outputs * (num_states + num_inputs) + num_states;
        for (auto i = 0u; i < num_states; ++i) {
            for (auto j = 0u; j < num_inputs; ++j) {
                B[i + j * num_states] = Jac0[count++];
            }
            count += num_states;

        }


        // return the struct holding all the results
        return linearization_results<value_type>{ std::move(A), std::move(B), std::move(C), std::move(D),
                                                  std::vector<value_type>(std::next(yxdot0.cbegin(), num_outputs), yxdot0.cend()),
                                                  std::vector<value_type>(yxdot0.cbegin(), std::next(yxdot0.cbegin(), num_outputs)),
                                                  std::vector<value_type>(xu0.cbegin(), std::next(xu0.cbegin(), num_states)),
                                                  std::vector<value_type>(std::next(xu0.cbegin(), num_states), xu0.cend()) };

    }


    template<typename ValueTypeOfResults = double,
             class Plant, typename StatesType, typename InputsType>
    inline linearization_results<ValueTypeOfResults> linearize_plant_numerically(const Plant &plant,
                                                                                 const StatesType &x0, const InputsType &u0,
                                                                                 std::vector<ValueTypeOfResults> deltas_numjac_states = { ValueTypeOfResults{ 1e-11 } },
                                                                                 std::vector<ValueTypeOfResults> deltas_numjac_inputs = { ValueTypeOfResults{ 1e-11 } })
    {
        //
        // Calculates the state-space representation
        // of a plant's linearization:
        //
        //  "y0   = plant.outputs(x0, u0),
        //  xdot0 = plant.derivatives(x0, y0, u0),
        //
        //  Dx    = x - x0,
        //  Du    = u - u0,
        //  Dy    = y - y0,
        //  Dxdot = xdot - xdot0,
        //
        //  Dxdot = A * Dx + B * Du,
        //  Dy    = C * Dx + D * Du",
        //
        // at a given point "[x0, u0]" with central
        // finite differences, and returns a
        // "linearization_results" object with
        // the results.
        //

        using value_type = ValueTypeOfResults;


        // get the sizes of the plant
        const auto num_states = plant.num_states();
        const auto num_outputs = plant.num_outputs();
        const auto num_inputs = plant.num_inputs();


        // ensure sizes
        ensure_minimum_size(num_states, x0);
        ensure_minimum_size(num_inputs, u0);

        if (deltas_numjac_states.size() == 1u) {
            deltas_numjac_states.insert(deltas_numjac_states.end(), num_states - 1u,
                                        deltas_numjac_states.front());

        }
        else if (deltas_numjac_states.size() != num_states) {
            TR7_THROW_RUNTIME_ERROR("F16_Nguyen::linearize_plant_numerically: invalid size of"
                                    " input \"deltas_numjac_states\".")

        }

        if (deltas_numjac_inputs.size() == 1u) {
            deltas_numjac_inputs.insert(deltas_numjac_inputs.end(), num_inputs - 1u,
                                        deltas_numjac_inputs.front());

        }
        else if (deltas_numjac_inputs.size() != num_inputs) {
            TR7_THROW_RUNTIME_ERROR("F16_Nguyen::linearize_plant_numerically: invalid size of"
                                    " input \"deltas_numjac_inputs\".")

        }


        // copy states and inputs (to return them afterwards)
        std::vector<value_type> x(x0.begin(), x0.end());
        std::vector<value_type> u(u0.begin(), u0.end());


        // state matrix "A = dxdot/dx"; output matrix "C = dy/dx"
        std::vector<value_type> A(num_states * num_states);
        std::vector<value_type> C(num_outputs * num_states);
        for (auto j = 0u; j < num_states; ++j) {
            const auto X = x[j];

            x[j] = X + deltas_numjac_states[j];
            const auto y_p = plant.outputs(x, u);
            const auto xdot_p = plant.derivatives(x, y_p, u);

            x[j] = X - deltas_numjac_states[j];
            const auto y_m = plant.outputs(x, u);
            const auto xdot_m = plant.derivatives(x, y_m, u);

            x[j] = X;

            const auto den = value_type{ 1 } / (value_type{ 2 } * deltas_numjac_states[j]);
            for (auto i = 0u; i < num_states; ++i) {
                A[i + j * num_states] = den * (xdot_p[i] - xdot_m[i]);
            }
            for (auto i = 0u; i < num_outputs; ++i) {
                C[i + j * num_outputs] = den * (y_p[i] - y_m[i]);
            }

        }


        // input matrix "B = dxdot/du"; feedthrough matrix "D = dy/du"
        std::vector<value_type> B(num_states * num_inputs);
        std::vector<value_type> D(num_outputs * num_inputs);
        for (auto j = 0u; j < num_inputs; ++j) {
            const auto U = u[j];

            u[j] = U + deltas_numjac_inputs[j];
            const auto y_p = plant.outputs(x, u);
            const auto xdot_p = plant.derivatives(x, y_p, u);

            u[j] = U - deltas_numjac_inputs[j];
            const auto y_m = plant.outputs(x, u);
            const auto xdot_m = plant.derivatives(x, y_m, u);

            u[j] = U;

            const auto den = value_type{ 1 } / (value_type{ 2 } * deltas_numjac_inputs[j]);
            for (auto i = 0u; i < num_states; ++i) {
                B[i + j * num_states] = den * (xdot_p[i] - xdot_m[i]);
            }
            for (auto i = 0u; i < num_outputs; ++i) {
                D[i + j * num_outputs] = den * (y_p[i] - y_m[i]);
            }

        }


        // outputs & statesdots ant the linearization point
        const auto y = plant.outputs(x, u);
        const auto xdot = plant.derivatives(x, y, u);


        // return the struct holding all the results
        return linearization_results<value_type>{ std::move(A), std::move(B), std::move(C), std::move(D),
                                                  std::vector<value_type>(xdot.cbegin(), xdot.cend()),
                                                  std::vector<value_type>(y.cbegin(), y.cend()),
                                                  std::move(x),
                                                  std::move(u) };

    }


}


#endif