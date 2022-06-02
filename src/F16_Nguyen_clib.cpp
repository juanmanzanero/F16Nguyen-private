#include <cstring>

#include "F16_Nguyen/F16_Nguyen_clib.h"
#include "F16_Nguyen/F16_Nguyen_plant.h"
#include "F16_Nguyen/F16_Nguyen_trim.h"
#include "F16_Nguyen/linearize_plant.h"


//
// Helper functions.
//

template<typename Cont, typename It>
bool container2range(const Cont &c, It range_begin, std::size_t size_range)
{
    //
    // Copies a container onto the given range. If the size
    // of the range is greater than the size of the container, the
    // surplus elements of the range will be set to
    // "std::numeric_limits<double>::quiet_NaN()". If the size of the
    // range is less than that of the container, only "size_range"
    // elements will be copied. The function returns true if the
    // container and the range have the same size, false otherwise.
    //

    // check sizes
    bool samesize{ true };
    auto max_size{ c.size() };
    if (max_size > size_range) {
        samesize = false;
        max_size = size_range;

    }
    else if (max_size < size_range) {
        samesize = false;
        std::fill(range_begin + max_size, range_begin + size_range, std::numeric_limits<double>::quiet_NaN());

    }

    // copy
    std::copy_n(c.cbegin(), max_size, range_begin);
    return samesize;

}


template<typename ContainerOfstdStrings>
bool stringcontainer2cstrarray(char **csarr, std::size_t size_csarr,
                               const ContainerOfstdStrings &c)
{
    //
    // Copies the contents of a container of std::strings "c"
    // onto the array of C-strings "csarr" (preallocated to
    // size "size_csarr"), reallocating all of its char* elements.
    // If "size_sarr" is greater than the size of the container,
    // the surplus elements of "csarr" will be set to the empty
    // string. If "size_csarr" is less than that of the container,
    // only "size_csarr" elements will be copied. The function returns
    // true if "c.size() == size_csarr", false otherwise.
    //

    // check sizes
    bool samesize{ true };
    auto max_size{ c.size() };
    if (max_size > size_csarr) {
        samesize = false;
        max_size = size_csarr;

    }
    else if (max_size < size_csarr) {
        samesize = false;
        std::for_each(csarr + max_size, csarr + size_csarr,
                      [](auto &c) { if (c) { delete[] c; } c = new char[1]; *c = '\0'; });

    }

    // copy
    for (auto count = 0u; count < max_size; ++count) {
        if (csarr[count]) {
            delete[] csarr[count];
        }
        csarr[count] = new char[c[count].length() + 1u];

        std::strcpy(csarr[count], c[count].c_str());

    }

    return samesize;

}


template<std::size_t RequiredSizeForNoError>
void buffer_size_check(F16N_error_code &err,
                       std::size_t actual_size)
{
    //
    // Sets the input "F16N_error_code" to a value
    // that represents the result of checking
    // a given size against the one we really require.
    //

    if (actual_size != RequiredSizeForNoError) {
        err = F16N_bad_size;

        if (actual_size < RequiredSizeForNoError) {
            err = F16N_invalid_size;
        }

    }

}


template<class Names2StrFun>
F16N_error_code plant_names(char **names, int num_names, Names2StrFun names2str_fun)
{
    //
    // Fills the input array of char* with the strings in
    // the container returned by input function "names2str".
    // Input "names" MUST BE PREALLOCATED ON ENTRY to
    // size "num_names". If
    // "num_names < sizeof_container_returned_by_names2str", then
    // the state array will only get filled up to "num_names". If
    // "num_names > sizeof_container_returned_by_names2str", its
    // surplus elements will be set to the empty string.
    //
    // NOTE: each char* element of the "names" array will get allocated
    // inside this function, I've tested that this doesn't cause memory
    // leaks, Matlab is capable of deallocating "names" without problem.
    //

    if (!names) {
        return F16N_invalid_ptr;
    }

    return stringcontainer2cstrarray(names, num_names,
                                     names2str_fun()) ? F16N_success : F16N_bad_size;

}


template<std::size_t RequiredNumPlantStates,
         std::size_t RequiredNumPlantOutputs,
         std::size_t RequiredNumPlantInputs,
         typename ResultsType,
         class PlantCallerFun,
         typename StatesType, typename OutputsType, typename InputsType>
F16N_error_code call_plant(ResultsType *results, int num_results,
                           PlantCallerFun plant_caller_fun,
                           F16_Nguyen_clib_plant cplant,
                           const StatesType *states, int num_states,
                           const OutputsType *outputs, int num_outputs,
                           const InputsType *inputs, int num_inputs)
{
    //
    // Applies the input "plant_caller_fun" to
    // the "F16_Nguyen::F16_Nguyen_plant" hidden
    // inside the "cplant" input handle. This
    // function should have the signature
    //
    // "container_of_results = plant_caller_fun(plant, x, y, u)."
    //
    // where "plant" is a reference to the
    // "F16_Nguyen::F16_Nguyen_plant", and the rest
    // are const buffers. The returned
    // container of results will be copied onto the
    // "results" buffer, which should come preallocated
    // on entry to "num_results" (it may differ from
    // the size of the "container_of_results").
    //

    if (!cplant) {
        return F16N_invalid_plant;

    }
    else if (!results ||
             (RequiredNumPlantStates != 0u && !states) ||
             (RequiredNumPlantOutputs != 0u && !outputs) ||
             (RequiredNumPlantInputs != 0u && !inputs)) {

        if (results) {
            std::fill_n(results, num_results, std::numeric_limits<ResultsType>::quiet_NaN());
        }

        return F16N_invalid_ptr;

    }


    // check input buffer sizes, return with error code if
    // any of them is shorter than required
    F16N_error_code err = F16N_success;

    if constexpr (RequiredNumPlantStates != 0u) {
        buffer_size_check<RequiredNumPlantStates>(err, num_states);
        if (err == F16N_invalid_size) {
            std::fill_n(results, num_results, std::numeric_limits<ResultsType>::quiet_NaN());
            return err;

        }

    }

    if constexpr (RequiredNumPlantOutputs != 0u) {
        buffer_size_check<RequiredNumPlantOutputs>(err, num_outputs);
        if (err == F16N_invalid_size) {
            std::fill_n(results, num_results, std::numeric_limits<ResultsType>::quiet_NaN());
            return err;

        }

    }

    if constexpr (RequiredNumPlantInputs != 0u) {
        buffer_size_check<RequiredNumPlantInputs>(err, num_inputs);
        if (err == F16N_invalid_size) {
            std::fill_n(results, num_results, std::numeric_limits<ResultsType>::quiet_NaN());
            return err;

        }

    }


    // do the thing
    const auto &plant = *static_cast<const F16_Nguyen::F16_Nguyen_plant *>(cplant);
    auto states_ = F16_Nguyen::make_const_static_contiguous_range<RequiredNumPlantStates>(states);
    auto outputs_ = F16_Nguyen::make_const_static_contiguous_range<RequiredNumPlantOutputs>(outputs);
    auto inputs_ = F16_Nguyen::make_const_static_contiguous_range<RequiredNumPlantInputs>(inputs);

    return container2range(plant_caller_fun(plant, states_, outputs_, inputs_),
                           results, num_results) ? err : F16N_bad_size;

}


template<class PlantLinearizerFun>
F16N_error_code linearize_plant_impl(double *A, int size_A,
                                     double *B, int size_B,
                                     double *C, int size_C,
                                     double *D, int size_D,
                                     double *statedots, int num_statedots,
                                     double *outputs, int num_outputs,
                                     double *states, int num_states,
                                     double *inputs, int num_inputs,
                                     PlantLinearizerFun plant_linearizer_fun,
                                     F16_Nguyen_clib_plant cplant,
                                     const double *linearization_states, int num_linearization_states,
                                     const double *linearization_inputs, int num_linearization_inputs)
{
    //
    // Calculates the state-space representation
    // of the plant's linearization:
    //
    //  "y0   = outputs(x0, u0),
    //  xdot0 = derivatives(x0, y0, u0),
    //
    //  Dx    = x - x0,
    //  Du    = u - u0,
    //  Dy    = y - y0,
    //  Dxdot = xdot - xdot0,
    //
    //  Dxdot = A * Dx + B * Du,
    //  Dy    = C * Dx + D * Du",
    //
    // at a given point "[x0, u0]" (represented
    // by input buffers "linearization_states"
    // and "linearization_inputs", which should
    // come allocated with sizes
    // "num_linearization_states" and
    // "num_linearization_inputs") by applying
    // input function "plant_linearizer_fun",
    // which takes the plant, "x0", "u0" and returns
    // an "F16_Nguyen::linearization_results" struct.
    // The results are copied onto buffers "A", "B",
    // "C", "D", "statedots" (== "xdot0"), "outputs"
    // (== "y0"), "states" (== "x0") and "inputs"
    // (== "u0"), which should come preallocated
    // on entry to the given sizes. The sizes of
    // these resulting buffers may differ from the
    // actual ones of the plant.
    //

    //  helper function to set the results to NaN
    const auto results_to_NaN = [&]()
    {
        if (A) {
            std::fill_n(A, size_A, std::numeric_limits<double>::quiet_NaN());
        }
        if (B) {
            std::fill_n(B, size_B, std::numeric_limits<double>::quiet_NaN());
        }
        if (C) {
            std::fill_n(C, size_C, std::numeric_limits<double>::quiet_NaN());
        }
        if (D) {
            std::fill_n(D, size_D, std::numeric_limits<double>::quiet_NaN());
        }
        if (statedots) {
            std::fill_n(statedots, num_statedots, std::numeric_limits<double>::quiet_NaN());
        }
        if (outputs) {
            std::fill_n(outputs, num_outputs, std::numeric_limits<double>::quiet_NaN());
        }
        if (states) {
            std::fill_n(states, num_states, std::numeric_limits<double>::quiet_NaN());
        }
        if (inputs) {
            std::fill_n(inputs, num_inputs, std::numeric_limits<double>::quiet_NaN());
        }

    };


    // check input pointers
    if (!cplant) {
        return F16N_invalid_plant;

    }
    else if (!A || !B || !C || !D ||
             !statedots || !outputs || !states || !inputs ||
             !linearization_states || !linearization_inputs) {

        results_to_NaN();

        return F16N_invalid_ptr;

    }


    // check input buffer sizes, return with error code if
    // any of them is shorter than required
    constexpr auto RequiredNumPlantStates = F16_Nguyen::F16_Nguyen_plant::_num_states;
    constexpr auto RequiredNumPlantInputs = F16_Nguyen::F16_Nguyen_plant::_num_inputs;

    F16N_error_code err = F16N_success;

    buffer_size_check<RequiredNumPlantStates>(err, num_linearization_states);
    if (err == F16N_invalid_size) {
        results_to_NaN();
        return err;

    }

    buffer_size_check<RequiredNumPlantInputs>(err, num_linearization_inputs);
    if (err == F16N_invalid_size) {
        results_to_NaN();
        return err;

    }


    // do the thing
    const auto &plant = *static_cast<const F16_Nguyen::F16_Nguyen_plant *>(cplant);
    auto linearization_states_ = F16_Nguyen::make_const_static_contiguous_range<RequiredNumPlantStates>(linearization_states);
    auto linearization_inputs_ = F16_Nguyen::make_const_static_contiguous_range<RequiredNumPlantInputs>(linearization_inputs);

    const auto linearization_results = plant_linearizer_fun(plant,
                                                            linearization_states_,
                                                            linearization_inputs_);


    // copy the linearization results onto the output buffers
    err = container2range(linearization_results.A, A, size_A) ? err : F16N_bad_size;
    err = container2range(linearization_results.B, B, size_B) ? err : F16N_bad_size;
    err = container2range(linearization_results.C, C, size_C) ? err : F16N_bad_size;
    err = container2range(linearization_results.D, D, size_D) ? err : F16N_bad_size;
    err = container2range(linearization_results.xdot, statedots, num_statedots) ? err : F16N_bad_size;
    err = container2range(linearization_results.y, outputs, num_outputs) ? err : F16N_bad_size;
    err = container2range(linearization_results.x, states, num_states) ? err : F16N_bad_size;
    err = container2range(linearization_results.u, inputs, num_inputs) ? err : F16N_bad_size;

    return err;

}


//
// "F16_Nguyen_clib" impl.
//

int F16_Nguyen_clib_num_plant_states()
{
    //
    // Returns the number of states of
    // class "F16_Nguyen::F16_Nguyen_plant".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_plant::_num_states);

}


int F16_Nguyen_clib_num_plant_outputs()
{
    //
    // Returns the number of outputs of
    // class "F16_Nguyen::F16_Nguyen_plant".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_plant::_num_outputs);

}


int F16_Nguyen_clib_num_plant_inputs()
{
    //
    // Returns the number of inputs of
    // class "F16_Nguyen::F16_Nguyen_plant".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_plant::_num_inputs);

}


F16N_error_code F16_Nguyen_clib_plant_state_names(char **names, int num_names)
{
    //
    // Fills the state array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_plant"'s
    // states.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_plant::state_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_output_names(char **names, int num_names)
{
    //
    // Fills the output array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_plant"'s
    // outputs.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_plant::output_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_input_names(char **names, int num_names)
{
    //
    // Fills the input array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_plant"'s
    // inputs.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_plant::input_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_default_inputs(double *default_inputs, int num_default_inputs)
{
    //
    // Calls the "default_inputs" member function of
    // class "F16_Nguyen::F16_Nguyen_plant", and copies
    // the results onto the "default_inputs" buffer",
    // preallocated on entry to "num_default_inputs",
    // which may differ from the actual number of plant
    // inputs).
    //

    return container2range(F16_Nguyen::F16_Nguyen_plant::default_inputs(),
                           default_inputs, num_default_inputs) ? F16N_success : F16N_bad_size;

}


F16N_error_code F16_Nguyen_clib_new_plant(F16_Nguyen_clib_plant *cplant,
                                          const char *path_aerodataset,
                                          const char *path_enginedataset)
{
    //
    // Allocates a new "F16_Nguyen::F16_Nguyen_plant",
    // constructed with the input ctor parameters (which
    // can be nullptr if you want to call the default ctor),
    // and establishes an opaque "cplant" handle to it.
    //

    // destroy the input plant if not nullptr
    if (*cplant != nullptr) {
        if (F16_Nguyen_clib_delete_plant(cplant) != F16N_success) {
            return F16N_plant_initialization_error;
        }

    }


    // do the thing
    F16_Nguyen::F16_Nguyen_plant *plant = nullptr;
    try {
        // allocate a new plant, use the default ctor
        // if any of the inputs is a nullptr
        plant = !path_aerodataset || !path_enginedataset ? new F16_Nguyen::F16_Nguyen_plant() :
                                                           new F16_Nguyen::F16_Nguyen_plant(path_aerodataset, path_enginedataset);

        // connect the plant to the opaque handle
        *cplant = static_cast<F16_Nguyen_clib_plant>(plant);

    }
    catch (...) {
        delete plant;
        return F16N_plant_initialization_error;

    }

    return F16N_success;

}


F16N_error_code F16_Nguyen_clib_delete_plant(F16_Nguyen_clib_plant *cplant)
{
    //
    // Finalizes the input "cplant", by deleting
    // the "F16_Nguyen::F16_Nguyen_plant" object
    // it hides.
    //

    if (!*cplant) {
        return F16N_invalid_plant;
    }


    // release all memory & reset the handle
    auto *plant = static_cast<F16_Nguyen::F16_Nguyen_plant *>(*cplant);
    delete plant;
    *cplant = nullptr;

    return F16N_success;

}


F16N_error_code F16_Nguyen_clib_plant_build(F16_Nguyen_clib_plant cplant,
                                            const char *path_aerodataset,
                                            const char *path_enginedataset)
{
    //
    // Calls the "build" member function of the input
    // "F16_Nguyen::F16_Nguyen_plant" (hidden inside
    // the "cplant" handle).
    //

    if (!cplant) {
        return F16N_invalid_plant;
    }
    else if (!path_aerodataset || !path_enginedataset) {
        return F16N_invalid_ptr;
    }


    // call the "build" member function, reset the plant if we fail
    auto &plant = *static_cast<F16_Nguyen::F16_Nguyen_plant *>(cplant);
    try {
        plant.build(path_aerodataset, path_enginedataset);

    }
    catch (...) {
        plant = F16_Nguyen::F16_Nguyen_plant{};
        return F16N_plant_initialization_error;

    }

    return F16N_success;

}


F16N_error_code F16_Nguyen_clib_plant_outputs(double *outputs, int num_outputs,
                                              F16_Nguyen_clib_plant cplant,
                                              const double *states, int num_states,
                                              const double *inputs, int num_inputs)
{
    //
    // Calls the "outputs" member function of the
    // "F16_Nguyen::F16_Nguyen_plant" hidden inside
    // the input "cplant" handle. The function is
    // called with the "states" and "inputs" arguments,
    // which should be allocated with the sizes
    // indicated by inputs "num_states/inputs". The
    // results are copied onto the "outputs" buffer,
    // which should be preallocated on entry to
    // size "num_outputs" (which may differ from
    // the actual number of plant outputs).
    //

    const char *nullptr_ = nullptr;
    return call_plant<F16_Nguyen::F16_Nguyen_plant::_num_states,
                      0u,
                      F16_Nguyen::F16_Nguyen_plant::_num_inputs>(outputs, num_outputs,
                                                                 [](const auto &plant, const auto &x, const auto &, const auto &u) { return plant.outputs(x, u); },
                                                                 cplant,
                                                                 states, num_states,
                                                                 nullptr_, 0,
                                                                 inputs, num_inputs);

}


F16N_error_code F16_Nguyen_clib_plant_derivatives(double *statedots, int num_statedots,
                                                  F16_Nguyen_clib_plant cplant,
                                                  const double *states, int num_states,
                                                  const double *outputs, int num_outputs,
                                                  const double *inputs, int num_inputs)
{
    //
    // Calls the "derivatives" member function of the
    // "F16_Nguyen::F16_Nguyen_plant" hidden inside
    // the input "cplant" handle. The function is
    // called with the "states", "outputs" and "inputs"
    // arguments, which should be allocated with the sizes
    // indicated by inputs "num_states/outputs/inputs". The
    // results are copied onto the "statedots" buffer,
    // which should be preallocated on entry to
    // size "num_statedots" (which may differ from
    // the actual number of plant states).
    //

    return call_plant<F16_Nguyen::F16_Nguyen_plant::_num_states,
                      F16_Nguyen::F16_Nguyen_plant::_num_outputs,
                      F16_Nguyen::F16_Nguyen_plant::_num_inputs>(statedots, num_statedots,
                                                                 [](const auto &plant, const auto &x, const auto &y, const auto &u) { return plant.derivatives(x, y, u); },
                                                                 cplant,
                                                                 states, num_states,
                                                                 outputs, num_outputs,
                                                                 inputs, num_inputs);

}


F16N_error_code F16_Nguyen_clib_plant_state_limiters(double *limstates, int num_limstates,
                                                     F16_Nguyen_clib_plant cplant,
                                                     const double *states, int num_states,
                                                     const double *outputs, int num_outputs,
                                                     const double *inputs, int num_inputs)
{
    //
    // Calls the "state_limiters" member function of the
    // "F16_Nguyen::F16_Nguyen_plant" hidden inside
    // the input "cplant" handle. The function is
    // called with the "states", "outputs" and "inputs"
    // arguments, which should be allocated with the sizes
    // indicated by inputs "num_states/outputs/inputs". The
    // results are copied onto the "limstates" buffer,
    // which should be preallocated on entry to
    // size "num_limstates" (which may differ from
    // the actual number of plant states).
    //

    return call_plant<F16_Nguyen::F16_Nguyen_plant::_num_states,
                      F16_Nguyen::F16_Nguyen_plant::_num_outputs,
                      F16_Nguyen::F16_Nguyen_plant::_num_inputs>(limstates, num_limstates,
                                                                 [](const auto &plant, const auto &x, const auto &y, const auto &u) { return plant.state_limiters(x, y, u); },
                                                                 cplant,
                                                                 states, num_states,
                                                                 outputs, num_outputs,
                                                                 inputs, num_inputs);

}


int F16_Nguyen_clib_num_trim_inflags()
{
    //
    // Returns the number of inflags of
    // class "F16_Nguyen::F16_Nguyen_trim".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_trim::_num_inflags);

}


int F16_Nguyen_clib_num_trim_inputs()
{
    //
    // Returns the number of inputs of
    // class "F16_Nguyen::F16_Nguyen_trim".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_trim::_num_inputs);

}


int F16_Nguyen_clib_num_trim_outputs()
{
    //
    // Returns the number of outputs of
    // class "F16_Nguyen::F16_Nguyen_trim".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_trim::_num_outputs);

}


F16N_error_code F16_Nguyen_clib_trim_inflag_names(char **names, int num_names)
{
    //
    // Fills the state array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_trim"'s
    // inflags.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_trim::inflag_names2str(); });

}


F16N_error_code F16_Nguyen_clib_trim_input_names(char **names, int num_names)
{
    //
    // Fills the state array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_trim"'s
    // inputs.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_trim::input_names2str(); });

}


F16N_error_code F16_Nguyen_clib_trim_output_names(char **names, int num_names)
{
    //
    // Fills the state array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_trim"'s
    // outputs.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_trim::output_names2str(); });

}


F16N_error_code F16_Nguyen_clib_trim_default_inflags(int *default_inflags, int num_default_inflags)
{
    //
    // Calls the "default_inflags" member function of
    // class "F16_Nguyen::F16_Nguyen_trim", and copies
    // the results onto the "default_inflags" buffer",
    // preallocated on entry to "num_default_inflags",
    // which may differ from the actual number of trim
    // inflags).
    //

    return container2range(F16_Nguyen::F16_Nguyen_trim::default_inflags(),
                           default_inflags, num_default_inflags) ? F16N_success : F16N_bad_size;

}


F16N_error_code F16_Nguyen_clib_trim_default_inputs(double *default_inputs, int num_default_inputs)
{
    //
    // Calls the "default_inputs" member function of
    // class "F16_Nguyen::F16_Nguyen_trim", and copies
    // the results onto the "default_inputs" buffer",
    // preallocated on entry to "num_default_inputs",
    // which may differ from the actual number of trim
    // inputs).
    //

    return container2range(F16_Nguyen::F16_Nguyen_trim::default_inputs(),
                           default_inputs, num_default_inputs) ? F16N_success : F16N_bad_size;

}


F16N_error_code F16_Nguyen_clib_trim_plant(int *trim_success,
                                           double *trim_outputs, int num_trim_outputs,
                                           F16_Nguyen_clib_plant cplant,
                                           const int *trim_inflags, int num_trim_inflags,
                                           const double *trim_inputs, int num_trim_inputs)
{
    //
    // Calls the "trim_plant" member function of the
    // "F16_Nguyen::F16_Nguyen_trim" class, with the
    // "F16_Nguyen::F16_Nguyen_plant" hidden inside
    // the input "cplant" handle as an input.
    //

    const auto trimmer_fun = [&trim_success](const auto &plant,
                                             const auto &,
                                             const auto &inflags, const auto &inputs)
    {
        const auto [success, trim_outputs] = F16_Nguyen::F16_Nguyen_trim::trim_plant(plant, inflags, inputs);
        *trim_success = static_cast<int>(success);
        return trim_outputs;

    };

    const char *nullptr_ = nullptr;
    return call_plant<0u,
                      F16_Nguyen::F16_Nguyen_trim::_num_inflags,
                      F16_Nguyen::F16_Nguyen_trim::_num_inputs>(trim_outputs, num_trim_outputs,
                                                                trimmer_fun,
                                                                cplant,
                                                                nullptr_, 0,
                                                                trim_inflags, num_trim_inflags,
                                                                trim_inputs, num_trim_inputs);

}


F16N_error_code F16_Nguyen_clib_linearize_plant(double *A, int size_A,
                                                double *B, int size_B,
                                                double *C, int size_C,
                                                double *D, int size_D,
                                                double *statedots, int num_statedots,
                                                double *outputs, int num_outputs,
                                                double *states, int num_states,
                                                double *inputs, int num_inputs,
                                                F16_Nguyen_clib_plant cplant,
                                                const double *linearization_states, int num_linearization_states,
                                                const double *linearization_inputs, int num_linearization_inputs)
{
    //
    // Linearizes the "F16_Nguyen::F16_Nguyen_plant"
    // (hidden inside input handle "cplant")
    // using autodiff, calling "linearize_plant_impl".
    //

    return linearize_plant_impl(A, size_A,
                                B, size_B,
                                C, size_C,
                                D, size_D,
                                statedots, num_statedots,
                                outputs, num_outputs,
                                states, num_states,
                                inputs, num_inputs,
                                [&](const auto &plant, const auto &linearization_states_, const auto &linearization_inputs_)
                                {
                                    return F16_Nguyen::linearize_plant(plant,
                                                                       linearization_states_,
                                                                       linearization_inputs_);

                                },
                                cplant,
                                linearization_states, num_linearization_states,
                                linearization_inputs, num_linearization_inputs);

}


F16N_error_code F16_Nguyen_clib_linearize_plant_numerically(double *A, int size_A,
                                                            double *B, int size_B,
                                                            double *C, int size_C,
                                                            double *D, int size_D,
                                                            double *statedots, int num_statedots,
                                                            double *outputs, int num_outputs,
                                                            double *states, int num_states,
                                                            double *inputs, int num_inputs,
                                                            F16_Nguyen_clib_plant cplant,
                                                            const double *linearization_states, int num_linearization_states,
                                                            const double *linearization_inputs, int num_linearization_inputs,
                                                            const double *deltas_numjac_states, int num_deltas_numjac_states,
                                                            const double *deltas_numjac_inputs, int num_deltas_numjac_inputs)
{
    //
    // Linearizes the "F16_Nguyen::F16_Nguyen_plant"
    // (hidden inside input handle "cplant")
    // numerically, calling "linearize_plant_impl".
    //

    //  helper function to set the results to NaN
    const auto results_to_NaN = [&]()
    {
        if (A) {
            std::fill_n(A, size_A, std::numeric_limits<double>::quiet_NaN());
        }
        if (B) {
            std::fill_n(B, size_B, std::numeric_limits<double>::quiet_NaN());
        }
        if (C) {
            std::fill_n(C, size_C, std::numeric_limits<double>::quiet_NaN());
        }
        if (D) {
            std::fill_n(D, size_D, std::numeric_limits<double>::quiet_NaN());
        }
        if (statedots) {
            std::fill_n(statedots, num_statedots, std::numeric_limits<double>::quiet_NaN());
        }
        if (outputs) {
            std::fill_n(outputs, num_outputs, std::numeric_limits<double>::quiet_NaN());
        }
        if (states) {
            std::fill_n(states, num_states, std::numeric_limits<double>::quiet_NaN());
        }
        if (inputs) {
            std::fill_n(inputs, num_inputs, std::numeric_limits<double>::quiet_NaN());
        }

    };


    // check sizes
    if ((deltas_numjac_states && !deltas_numjac_inputs) ||
        (!deltas_numjac_states && deltas_numjac_inputs)) {

        results_to_NaN();
        return F16N_invalid_ptr;

    }

    if ((num_deltas_numjac_inputs == 0 && num_deltas_numjac_states != 0) ||
        (num_deltas_numjac_inputs != 0 && num_deltas_numjac_states == 0)) {

        results_to_NaN();
        return F16N_invalid_size;

    }

    if ((num_deltas_numjac_states > 1 && num_deltas_numjac_states < F16_Nguyen::F16_Nguyen_plant::_num_states) ||
        (num_deltas_numjac_inputs > 1 && num_deltas_numjac_inputs < F16_Nguyen::F16_Nguyen_plant::_num_inputs)) {

        results_to_NaN();
        return F16N_invalid_size;

    }


    // do the thing
    return linearize_plant_impl(A, size_A,
                                B, size_B,
                                C, size_C,
                                D, size_D,
                                statedots, num_statedots,
                                outputs, num_outputs,
                                states, num_states,
                                inputs, num_inputs,
                                [&](const auto &plant, const auto &linearization_states_, const auto &linearization_inputs_)
                                {
                                    if ((!deltas_numjac_states && !deltas_numjac_inputs) ||
                                        (num_deltas_numjac_states == 0 && num_deltas_numjac_inputs == 0)) {

                                        return F16_Nguyen::linearize_plant_numerically(plant,
                                                                                       linearization_states_,
                                                                                       linearization_inputs_);

                                    }
                                    else {
                                        return F16_Nguyen::linearize_plant_numerically(plant,
                                                                                       linearization_states_,
                                                                                       linearization_inputs_,
                                                                                       std::vector<double>(deltas_numjac_states, deltas_numjac_states + num_deltas_numjac_states),
                                                                                       std::vector<double>(deltas_numjac_inputs, deltas_numjac_inputs + num_deltas_numjac_inputs));

                                    }

                                },
                                cplant,
                                linearization_states, num_linearization_states,
                                linearization_inputs, num_linearization_inputs);

}


int F16_Nguyen_clib_num_plant_aerodataset_inputs()
{
    //
    // Returns the number of inputs to the
    // "F16_Nguyen::F16_Nguyen_plant"'s aero-dataset.
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_aerodynamics::_num_dataset_inputs);

}


int F16_Nguyen_clib_num_plant_aerodataset_coefficients()
{
    //
    // Returns the number of coefficients in the
    // "F16_Nguyen::F16_Nguyen_plant"'s aero-dataset.
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_aerodynamics::_num_dataset_coefficients);

}


F16N_error_code F16_Nguyen_clib_plant_aerodataset_input_names(char **names, int num_names)
{
    //
    // Fills the state array of char* with the names
    // of the inputs to the "F16_Nguyen::F16_Nguyen_plant"
    // aero-dataset.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_aerodynamics::dataset_input_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_aerodataset_coefficient_names(char **names, int num_names)
{
    //
    // Returns the names of the coefficients in the
    // "F16_Nguyen::F16_Nguyen_plant" aero-dataset.
    // See functions "F16_Nguyen_clib_plant_state_names",
    // "F16_Nguyen_clib_plant_output_names" or
    // "F16_Nguyen_clib_plant_input_names", this one
    // behaves exactly the same.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_aerodynamics::dataset_coefficient_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_aerodataset_coefficients(double *coeffs, int num_coeffs,
                                                               F16_Nguyen_clib_plant cplant,
                                                               const double *inputs, int num_inputs)
{
    //
    // Evaluates the coefficients in the aero-dataset
    // of the "F16_Nguyen::F16_Nguyen_plant" hidden inside
    // the input "cplant" handle.
    //

    const char *nullptr_ = nullptr;
    return call_plant<0u,
                      0u,
                      F16_Nguyen::F16_Nguyen_aerodynamics::_num_dataset_inputs>(coeffs, num_coeffs,
                                                                                [](const auto &plant, const auto &, const auto &, const auto &u_dataset) { return plant.aerodynamics().dataset_coefficients(u_dataset); },
                                                                                cplant,
                                                                                nullptr_, 0,
                                                                                nullptr_, 0,
                                                                                inputs, num_inputs);

}


int F16_Nguyen_clib_num_plant_aeroforce_and_moment_inputs()
{
    //
    // Returns the number of inputs to the
    // "F16_Nguyen::F16_Nguyen_plant"'s calculation
    // of the total aerodynamic force and moment
    // coefficients.
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_aerodynamics::_num_force_and_moment_inputs);

}


int F16_Nguyen_clib_num_plant_aeroforce_and_moment_coefficients()
{
    //
    // Returns the number of total aero-force & moment
    // coefficients in the "F16_Nguyen::F16_Nguyen_plant"'s
    // aerodynamic model.
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_aerodynamics::_num_force_and_moment_coefficients);

}


F16N_error_code F16_Nguyen_clib_plant_aeroforce_and_moment_input_names(char **names, int num_names)
{
    //
    // Fills the state array of char* with the names
    // of the inputs to the "F16_Nguyen::F16_Nguyen_plant"
    // calculation of total aero-force and moment coefficients.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_aerodynamics::force_and_moment_input_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_aeroforce_and_moment_coefficient_names(char **names, int num_names)
{
    //
    // Returns the names of the total aero-force &
    // moment coefficients that the aerodynamic model
    // of "F16_Nguyen::F16_Nguyen_plant" calculates.
    // See functions "F16_Nguyen_clib_plant_state_names",
    // "F16_Nguyen_clib_plant_output_names" or
    // "F16_Nguyen_clib_plant_input_names", this one
    // behaves exactly the same.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_aerodynamics::force_and_moment_coefficient_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_aeroforce_and_moment_coefficients(double *coeffs, int num_coeffs,
                                                                        F16_Nguyen_clib_plant cplant,
                                                                        const double *inputs, int num_inputs)
{
    //
    // Evaluates the total aero-force & moment coefficients
    // for the "F16_Nguyen::F16_Nguyen_plant" hidden inside
    // the input "cplant" handle.
    //

    const char *nullptr_ = nullptr;
    return call_plant<0u,
                      0u,
                      F16_Nguyen::F16_Nguyen_aerodynamics::_num_force_and_moment_inputs>(
                          coeffs, num_coeffs,
                          [](const auto &plant, const auto &, const auto &, const auto &u_FM_short) { return plant.aerodynamics().force_and_moment_coefficients(u_FM_short); },
                          cplant,
                          nullptr_, 0,
                          nullptr_, 0,
                          inputs, num_inputs);

}


int F16_Nguyen_clib_num_plant_engine_states()
{
    //
    // Returns the number of states of
    // class "F16_Nguyen::F16_Nguyen_engine".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_engine::_num_states);

}


int F16_Nguyen_clib_num_plant_engine_outputs()
{
    //
    // Returns the number of outputs of
    // class "F16_Nguyen::F16_Nguyen_engine".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_engine::_num_outputs);

}


int F16_Nguyen_clib_num_plant_engine_dataset_outputs()
{
    //
    // Returns the number of "dataset_outputs" of
    // class "F16_Nguyen::F16_Nguyen_engine".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_engine::_num_dataset_outputs);

}


int F16_Nguyen_clib_num_plant_engine_inputs()
{
    //
    // Returns the number of inputs of
    // class "F16_Nguyen::F16_Nguyen_engine".
    //

    return static_cast<int>(F16_Nguyen::F16_Nguyen_engine::_num_inputs);

}


F16N_error_code F16_Nguyen_clib_plant_engine_state_names(char **names, int num_names)
{
    //
    // Fills the state array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_engine"'s
    // states.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_engine::state_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_engine_output_names(char **names, int num_names)
{
    //
    // Fills the output array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_engine"'s
    // outputs.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_engine::output_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_engine_dataset_output_names(char **names, int num_names)
{
    //
    // Fills the output array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_engine"'s
    // "dataset_outputs".
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_engine::dataset_output_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_engine_input_names(char **names, int num_names)
{
    //
    // Fills the input array of char* with the names of all
    // of the scalar signals of "F16_Nguyen::F16_Nguyen_engine"'s
    // inputs.
    //

    return plant_names(names, num_names,
                       []() { return F16_Nguyen::F16_Nguyen_engine::input_names2str(); });

}


F16N_error_code F16_Nguyen_clib_plant_engine_outputs(double *outputs, int num_outputs,
                                                     F16_Nguyen_clib_plant cplant,
                                                     const double *states, int num_states,
                                                     const double *inputs, int num_inputs)
{
    //
    // Calls the "outputs" member function of the
    // "F16_Nguyen::F16_Nguyen_engine" hidden inside
    // the input "cplant" handle. The function is
    // called with the "states" and "inputs" arguments,
    // which should be allocated with the sizes
    // indicated by inputs "num_states/inputs". The
    // results are copied onto the "outputs" buffer,
    // which should be preallocated on entry to
    // size "num_outputs" (which may differ from
    // the actual number of plant outputs).
    //

    const char *nullptr_ = nullptr;
    return call_plant<F16_Nguyen::F16_Nguyen_engine::_num_states,
                      0u,
                      F16_Nguyen::F16_Nguyen_engine::_num_inputs>(outputs, num_outputs,
                                                                 [](const auto &plant, const auto &x, const auto &, const auto &u) { return plant.engine().outputs(x, u); },
                                                                 cplant,
                                                                 states, num_states,
                                                                 nullptr_, 0,
                                                                 inputs, num_inputs);

}


F16N_error_code F16_Nguyen_clib_plant_engine_dataset(double *dataset_outputs, int num_dataset_outputs,
                                                     F16_Nguyen_clib_plant cplant,
                                                     const double *inputs, int num_inputs)
{
    //
    // Calls the "dataset" member function of the
    // "F16_Nguyen::F16_Nguyen_engine" hidden inside
    // the input "cplant" handle. The function is
    // called with the "inputs" arguments,
    // which should be allocated with the sizes
    // indicated by argument "num_inputs". The
    // results are copied onto the "dataset_outputs"
    // buffer, which should be preallocated on entry to
    // size "num_dataset_outputs" (which may differ from
    // the actual number of plant outputs).
    //

    const char *nullptr_ = nullptr;
    return call_plant<0u,
                      0u,
                      F16_Nguyen::F16_Nguyen_engine::_num_inputs>(dataset_outputs, num_dataset_outputs,
                                                                 [](const auto &plant, const auto &, const auto &, const auto &u) { return plant.engine().dataset(u); },
                                                                 cplant,
                                                                 nullptr_, 0,
                                                                 nullptr_, 0,
                                                                 inputs, num_inputs);

}


F16N_error_code F16_Nguyen_clib_plant_engine_derivatives(double *statedots, int num_statedots,
                                                         F16_Nguyen_clib_plant cplant,
                                                         const double *states, int num_states,
                                                         const double *outputs, int num_outputs)
{
    //
    // Calls the "derivatives" member function of the
    // "F16_Nguyen::F16_Nguyen_engine" hidden inside
    // the input "cplant" handle. The function is
    // called with the "states", "outputs" and "inputs"
    // arguments, which should be allocated with the sizes
    // indicated by inputs "num_states/outputs/inputs". The
    // results are copied onto the "statedots" buffer,
    // which should be preallocated on entry to
    // size "num_statedots" (which may differ from
    // the actual number of plant states).
    //

    const char *nullptr_ = nullptr;
    return call_plant<F16_Nguyen::F16_Nguyen_engine::_num_states,
                      F16_Nguyen::F16_Nguyen_engine::_num_outputs,
                      0u>(statedots, num_statedots,
                          [](const auto &plant, const auto &x, const auto &y, const auto &) { return plant.engine().derivatives(x, y); },
                          cplant,
                          states, num_states,
                          outputs, num_outputs,
                          nullptr_, 0);

}


F16N_error_code F16_Nguyen_clib_plant_engine_state_limiters(double *limstates, int num_limstates,
                                                            F16_Nguyen_clib_plant cplant,
                                                            const double *states, int num_states)
{
    //
    // Calls the "state_limiters" member function of the
    // "F16_Nguyen::F16_Nguyen_engine" hidden inside
    // the input "cplant" handle. The function is
    // called with the "states", "outputs" and "inputs"
    // arguments, which should be allocated with the sizes
    // indicated by inputs "num_states/outputs/inputs". The
    // results are copied onto the "limstates" buffer,
    // which should be preallocated on entry to
    // size "num_limstates" (which may differ from
    // the actual number of plant states).
    //

    const char *nullptr_ = nullptr;
    return call_plant<F16_Nguyen::F16_Nguyen_engine::_num_states,
                      0u,
                      0u>(limstates, num_limstates,
                          [](const auto &plant, const auto &x, const auto &, const auto &) { return plant.engine().state_limiters(x); },
                          cplant,
                          states, num_states,
                          nullptr_, 0,
                          nullptr_, 0);

}
