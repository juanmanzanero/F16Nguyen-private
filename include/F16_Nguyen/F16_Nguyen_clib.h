#ifndef F16_NGUYEN_CLIB_H
#define F16_NGUYEN_CLIB_H


//
// A C interface to call "F16_Nguyen".
//


//
// Windows/Linux API declaration (exported functions).
//

#ifdef _MSC_VER
#ifdef F16_NGUYEN_CLIB_EXPORTS
#define F16_NGUYEN_CLIB_API __declspec(dllexport)
#else
#define F16_NGUYEN_CLIB_API __declspec(dllimport)
#endif

#elif (defined (__GNUC__) || defined (__GNUG__))
#ifdef F16_NGUYEN_CLIB_EXPORTS
// remember to compile with -fvisibility=hidden
#define F16_NGUYEN_CLIB_API __attribute__((visibility ("default")))
#else
#define F16_NGUYEN_CLIB_API
#endif

#else
#define F16_NGUYEN_CLIB_API
#endif


//
// Error handling.
//

typedef enum
{
    F16N_success = 0,

    // irrecoverable errors (< 0), something's clearly wrong
    F16N_plant_initialization_error = -1,
    F16N_invalid_plant = -2,

    // recoverable errors (> 0), try again changing your call
    F16N_invalid_ptr = 3,
    F16N_invalid_size = 3,
    F16N_bad_size = 4,

} F16N_errors;

typedef int F16N_error_code;


//
// Plant instance: an opaque handle
// that a user of the library will own,
// so that the library can have many users,
// each of them calling it with their own
// instance.
//

typedef void *F16_Nguyen_clib_plant;


//
// Library API. We also provide a fun-pointer typedef
// per function, so that users that open the library
// with dlopen can write cleaner code.
//

#ifdef __cplusplus
extern "C"
{
#endif

    //
    // Query the number & names of the states,
    // outputs and inputs of the plant.
    //

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_states();
    typedef int(*F16_Nguyen_clib_num_plant_states_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_outputs();
    typedef int(*F16_Nguyen_clib_num_plant_outputs_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_inputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_plant_inputs_handle)();

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_state_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_state_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_output_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_output_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_input_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_input_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_default_inputs(double *default_inputs, int num_default_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_default_inputs_handle)(double *, int);


    //
    // Create/rebuild/destroy a plant.
    //

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_new_plant(F16_Nguyen_clib_plant *cplant,
                                                  const char *path_aerodataset,
                                                  const char *path_enginedataset);
    typedef F16N_error_code(*F16_Nguyen_clib_new_plant_handle)(F16_Nguyen_clib_plant *,
                                                               const char *,
                                                               const char *);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_delete_plant(F16_Nguyen_clib_plant *cplant);
    typedef F16N_error_code(*F16_Nguyen_clib_delete_plant_handle)(F16_Nguyen_clib_plant *);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_build(F16_Nguyen_clib_plant cplant,
                                                    const char *path_aerodataset,
                                                    const char *path_enginedataset);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_build_handle)(F16_Nguyen_clib_plant,
                                                                 const char *,
                                                                 const char *);


    //
    // Call the plant's state-space functions.
    //

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_outputs(double *outputs, int num_outputs,
                                                      F16_Nguyen_clib_plant cplant,
                                                      const double *states, int num_states,
                                                      const double *inputs, int num_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_outputs_handle)(double *, int,
                                                                   F16_Nguyen_clib_plant,
                                                                   const double *, int,
                                                                   const double *, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_derivatives(double *statedots, int num_statedots,
                                                          F16_Nguyen_clib_plant cplant,
                                                          const double *states, int num_states,
                                                          const double *outputs, int num_outputs,
                                                          const double *inputs, int num_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_derivatives_handle)(double *, int,
                                                                       F16_Nguyen_clib_plant,
                                                                       const double *, int,
                                                                       const double *, int,
                                                                       const double *, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_state_limiters(double *stateslim, int num_stateslim,
                                                             F16_Nguyen_clib_plant cplant,
                                                             const double *states, int num_states,
                                                             const double *outputs, int num_outputs,
                                                             const double *inputs, int num_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_state_limiters_handle)(double *, int,
                                                                          F16_Nguyen_clib_plant,
                                                                          const double *, int,
                                                                          const double *, int,
                                                                          const double *, int);


    //
    // Trim the plant.
    //

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_trim_inflags();
    typedef F16N_error_code(*F16_Nguyen_clib_num_trim_inflags_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_trim_inputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_trim_inputs_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_trim_outputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_trim_outputs_handle)();

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_trim_inflag_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_trim_inflag_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_trim_input_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_trim_input_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_trim_output_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_trim_output_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_trim_default_inflags(int *default_inflags, int num_default_inflags);
    typedef F16N_error_code(*F16_Nguyen_clib_trim_default_inflags_handle)(int *, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_trim_default_inputs(double *default_inputs, int num_default_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_trim_default_inputs_handle)(double *, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_trim_plant(int *trim_success,
                                                   double *trim_outputs, int num_trim_outputs,
                                                   F16_Nguyen_clib_plant cplant,
                                                   const int *trim_inflags, int num_trim_inflags,
                                                   const double *trim_inputs, int num_trim_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_trim_plant_handle)(int *,
                                                                double *, int,
                                                                F16_Nguyen_clib_plant,
                                                                const int*, int,
                                                                const double *, int);


    //
    // Linearize the plant.
    //

    extern F16_NGUYEN_CLIB_API
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
                                                        const double *linearization_inputs, int num_linearization_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_linearize_plant_handle)(double *, int,
                                                                     double *, int,
                                                                     double *, int,
                                                                     double *, int,
                                                                     double *, int,
                                                                     double *, int,
                                                                     double *, int,
                                                                     double *, int,
                                                                     F16_Nguyen_clib_plant,
                                                                     const double *, int,
                                                                     const double *, int);

    extern F16_NGUYEN_CLIB_API
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
                                                                    const double *deltas_numjac_inputs, int num_deltas_numjac_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_linearize_plant_numerically_handle)(double *, int,
                                                                                 double *, int,
                                                                                 double *, int,
                                                                                 double *, int,
                                                                                 double *, int,
                                                                                 double *, int,
                                                                                 double *, int,
                                                                                 double *, int,
                                                                                 F16_Nguyen_clib_plant,
                                                                                 const double *, int,
                                                                                 const double *, int,
                                                                                 const double *, int,
                                                                                 const double *, int);


    //
    // Call the plant's aero-dataset (i.e., the lookup tables
    // that provide all of the aero model's coefficients).
    //

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_aerodataset_inputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_aerodataset_inputs_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_aerodataset_coefficients();
    typedef F16N_error_code(*F16_Nguyen_clib_num_aerodataset_coefficients_handle)();

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_aerodataset_input_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_aerodataset_input_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_aerodataset_coefficient_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_aerodataset_coefficient_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_aerodataset_coefficients(double *coeffs, int num_coeffs,
                                                                       F16_Nguyen_clib_plant cplant,
                                                                       const double *inputs, int num_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_aerodataset_coefficients_handle)(double *, int,
                                                                              F16_Nguyen_clib_plant,
                                                                              const double *, int);


    //
    // Call the plant's aerodynamic model, which yields
    // the total aero-force and moment coefficients.
    //

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_aeroforce_and_moment_inputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_plant_aeroforce_and_moment_inputs_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_aeroforce_and_moment_coefficients();
    typedef F16N_error_code(*F16_Nguyen_clib_num_plant_aeroforce_and_moment_coefficients_handle)();

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_aeroforce_and_moment_input_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_aeroforce_and_moment_input_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_aeroforce_and_moment_coefficient_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_aeroforce_and_moment_coefficient_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_aeroforce_and_moment_coefficients(double *coeffs, int num_coeffs,
                                                                                F16_Nguyen_clib_plant cplant,
                                                                                const double *inputs, int num_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_aeroforce_and_moment_coefficients_handle)(double *, int,
                                                                                       F16_Nguyen_clib_plant,
                                                                                       const double *, int);


    //
    // Call the plant's engine model, which yields
    // the thrust and the engine state.
    //

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_engine_states();
    typedef F16N_error_code(*F16_Nguyen_clib_num_num_engine_states_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_engine_outputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_num_engine_outputs_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_engine_dataset_outputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_num_engine_dataset_outputs_handle)();

    extern F16_NGUYEN_CLIB_API
        int F16_Nguyen_clib_num_plant_engine_inputs();
    typedef F16N_error_code(*F16_Nguyen_clib_num_num_engine_inputs_handle)();

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_state_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_state_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_output_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_output_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_input_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_input_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_dataset_output_names(char **names, int num_names);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_dataset_output_names_handle)(char **, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_outputs(double *outputs, int num_outputs,
                                                             F16_Nguyen_clib_plant cplant,
                                                             const double *states, int num_states,
                                                             const double *inputs, int num_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_outputs_handle)(double *, int,
                                                                          F16_Nguyen_clib_plant,
                                                                          const double *, int,
                                                                          const double *, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_dataset(double *dataset_outputs, int num_dataset_outputs,
                                                             F16_Nguyen_clib_plant cplant,
                                                             const double *inputs, int num_inputs);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_dataset_handle)(double *, int,
                                                                          F16_Nguyen_clib_plant,
                                                                          const double *, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_derivatives(double *statedots, int num_statedots,
                                                                 F16_Nguyen_clib_plant cplant,
                                                                 const double *states, int num_states,
                                                                 const double *outputs, int num_outputs);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_derivatives_handle)(double *, int,
                                                                              F16_Nguyen_clib_plant,
                                                                              const double *, int,
                                                                              const double *, int);

    extern F16_NGUYEN_CLIB_API
        F16N_error_code F16_Nguyen_clib_plant_engine_state_limiters(double *stateslim, int num_stateslim,
                                                                    F16_Nguyen_clib_plant cplant,
                                                                    const double *states, int num_states);
    typedef F16N_error_code(*F16_Nguyen_clib_plant_engine_state_limiters_handle)(double *, int,
                                                                                 F16_Nguyen_clib_plant,
                                                                                 const double *, int);


#ifdef __cplusplus
}
#endif


#endif
