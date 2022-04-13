#ifndef F16_NGUYEN_F16_NGUYEN_AERODYNAMICS_H
#define F16_NGUYEN_F16_NGUYEN_AERODYNAMICS_H
#pragma once


#include "tr7/lookup.h"
#include "tr7/tr7math.h"

#include "F16_Nguyen/modelglue.h"
#include "F16_Nguyen/contiguous_ranges.h"
#include "F16_Nguyen/read_csv_table.h"


//
// Defines the "F16_Nguyen_aerodynamics" class, which implements
// the aerodynamic model explained in reference "Simulator study
// of stall/post-stall characteristics of a fighter airplane with
// relaxed longitudinal static stability", L.T.Nguyen et al.,
// NASA-TP-1538, 1979.
//


namespace F16_Nguyen
{
    class F16_Nguyen_aerodynamics
    {
        //
        // A class to calculate the aerodynamic
        // forces & moments of Nguyen's F16
        // simulation model.
        //

    public:

        using data_type = double;


        struct default_parameters
        {
            static constexpr data_type MAC_m = tr7::ft2m(11.32);
            static constexpr data_type wingspan_m = tr7::ft2m(30.);
            static constexpr data_type wing_surface_m2 = tr7::ft2m(tr7::ft2m(300.));

            static constexpr data_type reference_xcg_per_MAC = 0.35;
            static constexpr data_type reference_ycg_per_semiwingspan = 0.;
            static constexpr data_type reference_leading_edge_flap_deg = 25.;
            static constexpr data_type reference_speedbreak_deg = 60.;
            static constexpr data_type reference_aileron_deg = 20.;
            static constexpr data_type reference_rudder_deg = 30.;

        };


        DECLARE_MODEL_SIGNALS(force_and_moment_coefficient, CXtot, CYtot, CZtot,
                                                            CLLtot, CMtot, CNtot)

        DECLARE_MODEL_SIGNALS(dataset_input, aoa_deg, aos_deg, dh_deg)

        DECLARE_MODEL_SIGNALS_DERIVED_FROM(force_and_moment_input, dataset_input, TAS_mps,
                                                                                  paero_radps, qaero_radps, raero_radps,
                                                                                  dlef_deg, dsb_deg,
                                                                                  da_deg, dr_deg,
                                                                                  aoadot_radps,
                                                                                  xcg_per_MAC, ycg_per_semiwingspan, MAC_m, wingspan_m)

        DECLARE_MODEL_OUTPUTS_DERIVED_FROM(force_and_moment_coefficient, Faerox_N, Faeroy_N, Faeroz_N,
                                                                         Maerox_Nm, Maeroy_Nm, Maeroz_Nm)

        DECLARE_MODEL_SIGNALS(dataset_coefficient, CX, CXdh0, CXlef, CXq, DCXsb, DCXqlef,
                                                   CY, CYda20, CYdr30, CYlef, CYda20lef, CYp, CYr, DCYplef, DCYrlef,
                                                   CZ, CZdh0, CZlef, CZq, DCZsb, DCZqlef,
                                                   CLL, CLLdh0, CLLda20, CLLdr30, CLLlef, CLLda20lef, CLLp, CLLr, DCLLbeta, DCLLplef, DCLLrlef,
                                                   CM, CMdh0, CMlef, DCMds, CMq, DCM, DCMsb, DCMqlef,
                                                   CN, CNdh0, CNda20, CNdr30, CNlef, CNda20lef, CNp, CNr, DCNbeta, DCNda, DCNplef, DCNrlef,
                                                   etadh,
                                                   CLalphadot, CMalphadot)


        F16_Nguyen_aerodynamics() = default;

        template<typename... Args>
        F16_Nguyen_aerodynamics(Args&&... args)
        {
            build(std::forward<Args>(args)...);
        }

        void build(const std::string &path_dataset);


        template<typename T>
        constexpr outputs_type<T> outputs(const T &dynamic_pressure_Pa,
                                          const T &TAS_mps, const T &aoa_deg, const T &aos_deg,
                                          const T &paero_radps, const T &qaero_radps, const T &raero_radps,
                                          const T &dh_deg, const T &dlef_deg, const T &dsb_deg,
                                          const T &da_deg, const T &dr_deg,
                                          const T &aoadot_radps,
                                          const T &xcg_per_MAC = T{ default_parameters::reference_xcg_per_MAC },
                                          const T &ycg_per_semiwingspan = T{ default_parameters::reference_ycg_per_semiwingspan },
                                          const T &MAC_m = T{ default_parameters::MAC_m },
                                          const T &wingspan_m = T{ default_parameters::wingspan_m },
                                          const T &wing_surface_m2 = T{ default_parameters::wing_surface_m2 }) const;


        template<typename T>
        constexpr force_and_moment_coefficients_type<T> force_and_moment_coefficients(const T &aoa_deg, const T &aos_deg,
                                                                                      const T &dh_deg,
                                                                                      const T &TAS_mps,
                                                                                      const T &paero_radps, const T &qaero_radps, const T &raero_radps,
                                                                                      const T &dlef_deg, const T &dsb_deg,
                                                                                      const T &da_deg, const T &dr_deg,
                                                                                      const T &aoadot_radps,
                                                                                      const T &xcg_per_MAC = T{ default_parameters::reference_xcg_per_MAC },
                                                                                      const T &ycg_per_semiwingspan = T{ default_parameters::reference_ycg_per_semiwingspan },
                                                                                      const T &MAC_m = T{ default_parameters::MAC_m },
                                                                                      const T &wingspan_m = T{ default_parameters::wingspan_m }) const;

        template<typename ForceAndMomentInputsType>
        constexpr auto force_and_moment_coefficients(const ForceAndMomentInputsType & u_FM) const
        {
            ensure_minimum_size<_num_force_and_moment_inputs>(u_FM);
            return force_and_moment_coefficients(u_FM[force_and_moment_input_names::aoa_deg], u_FM[force_and_moment_input_names::aos_deg],
                                                 u_FM[force_and_moment_input_names::dh_deg],
                                                 u_FM[force_and_moment_input_names::TAS_mps],
                                                 u_FM[force_and_moment_input_names::paero_radps], u_FM[force_and_moment_input_names::qaero_radps], u_FM[force_and_moment_input_names::raero_radps],
                                                 u_FM[force_and_moment_input_names::dlef_deg], u_FM[force_and_moment_input_names::dsb_deg],
                                                 u_FM[force_and_moment_input_names::da_deg], u_FM[force_and_moment_input_names::dr_deg],
                                                 u_FM[force_and_moment_input_names::aoadot_radps],
                                                 u_FM[force_and_moment_input_names::xcg_per_MAC],
                                                 u_FM[force_and_moment_input_names::ycg_per_semiwingspan],
                                                 u_FM[force_and_moment_input_names::MAC_m],
                                                 u_FM[force_and_moment_input_names::wingspan_m]);

        }


        template<typename T>
        constexpr dataset_coefficients_type<T> dataset_coefficients(const T &aoa_deg, const T &aos_deg,
                                                                    const T &dh_deg) const;

        template<typename DatasetInputsType>
        constexpr auto dataset_coefficients(const DatasetInputsType &u_dataset) const
        {
            ensure_minimum_size<_num_dataset_inputs>(u_dataset);
            return dataset_coefficients(u_dataset[dataset_input_names::aoa_deg], u_dataset[dataset_input_names::aos_deg],
                                        u_dataset[dataset_input_names::dh_deg]);

        }


    private:

        template<typename T>
        constexpr dataset_coefficients_type<T> lookup_dataset_coefficients(const T &aoa_deg, const T &aos_deg,
                                                                           const T &dh_deg) const;


        struct prelookups
        {
            static constexpr auto prelookup_option = tr7::lookup::prelookup::linsearch;
            static constexpr auto use_previous_index = true;

            using prelookup_type = tr7::Prelookup<prelookup_option,
                                                  use_previous_index,
                                                  std::vector<data_type> >;

            prelookup_type aoa_long_deg;
            prelookup_type aoa_short_deg;
            prelookup_type aos_deg;
            prelookup_type dh_long_deg;
            prelookup_type dh_mid_deg;
            prelookup_type dh_short_deg;

        } _prelookups;


        struct direct_luts
        {
            static constexpr auto interpolation_option = tr7::lookup::interpolation::linear;
            static constexpr auto extrapolation_option = tr7::lookup::extrapolation::nearest;

            template<std::size_t Rank, std::size_t Sizeof_Value>
            using direct_lut_type = tr7::Direct_lookup_table<Rank,
                                                             interpolation_option,
                                                             extrapolation_option,
                                                             std::vector<std::array<data_type, Sizeof_Value> > >;

            // coefficients that are "f(aoa_long, aos, dh_mid)"
            struct C_aoa_long_aos_dh_mid : direct_lut_type<3u, 3u>
            {
                enum { CX, CZ, CM };
            } C_aoa_long_aos_dh_mid;

            // coefficients that are "f(aoa_long, aos, dh_short)"
            struct C_aoa_long_aos_dh_short : direct_lut_type<3u, 2u>
            {
                enum { CLL, CN };
            } C_aoa_long_aos_dh_short;

            // coefficients that are "f(aoa_long, aos)"
            struct C_aoa_long_aos : direct_lut_type<2u, 12u>
            {
                enum { CXdh0, CY, CYda20, CYdr30, CZdh0, CLLdh0, CLLda20, CLLdr30, CMdh0, CNdh0, CNda20, CNdr30 };
            } C_aoa_long_aos;

            // coefficients that are "f(aoa_short, aos)"
            struct C_aoa_short_aos : direct_lut_type<2u, 9u>
            {
                enum { CXlef, CYlef, CYda20lef, CZlef, CLLlef, CLLda20lef, CMlef, CNlef, CNda20lef };
            } C_aoa_short_aos;

            // coefficients that are "f(aoa_long, dh_long)"
            struct C_aoa_long_dh_long : direct_lut_type<2u, 1u>
            {
                enum { DCMds };
            } C_aoa_long_dh_long;

            // coefficients that are "f(aoa_long)"
            struct C_aoa_long : direct_lut_type<1u, 16u>
            {
                enum { CXq, DCXsb, CYp, CYr, CZq, DCZsb, CLLp, CLLr, DCLLbeta, CMq, DCM, DCMsb, CNp, CNr, DCNbeta, DCNda };
            } C_aoa_long;

            // coefficients that are "f(aoa_short)"
            struct C_aoa_short : direct_lut_type<1u, 11u>
            {
                enum { DCXqlef, DCYplef, DCYrlef, DCZqlef, DCLLplef, DCLLrlef, DCMqlef, DCNplef, DCNrlef, CLalphadot, CMalphadot };
            } C_aoa_short;

            // coefficients that are "f(dh_mid)"
            struct C_dh_mid : direct_lut_type<1u, 1u>
            {
                enum { etadh };
            } C_dh_mid;

        } _dluts;

    };


    //
    // "F16_Nguyen_aerodynamics" impl.
    //

    inline void F16_Nguyen_aerodynamics::build(const std::string &path_dataset)
    {
        //
        // Reads the "NGUYEN F16 AERO-DATASET" "*.csv" files living in
        // "path_dataset", and sets up the lookup tables that we'll use
        // to evaluate the aircraft's aerodynamic coefficients.
        //

        // read breakpoints
        _prelookups.aoa_long_deg.grid_vector(std::get<0>(read_csv_table(path_dataset + "/aoa_long_breaks_deg.csv")));
        _prelookups.aoa_short_deg.grid_vector(std::get<0>(read_csv_table(path_dataset + "/aoa_short_breaks_deg.csv")));
        _prelookups.aos_deg.grid_vector(std::get<0>(read_csv_table(path_dataset + "/aos_breaks_deg.csv")));
        _prelookups.dh_long_deg.grid_vector(std::get<0>(read_csv_table(path_dataset + "/dh_long_breaks_deg.csv")));
        _prelookups.dh_mid_deg.grid_vector(std::get<0>(read_csv_table(path_dataset + "/dh_mid_breaks_deg.csv")));
        _prelookups.dh_short_deg.grid_vector(std::get<0>(read_csv_table(path_dataset + "/dh_short_breaks_deg.csv")));


        // read tables for "C_aoa_long_aos_dh_mid"
        {
            const auto size_x = _prelookups.aoa_long_deg.grid_vector().size();
            const auto size_y = _prelookups.aos_deg.grid_vector().size();
            const auto size_z = _prelookups.dh_mid_deg.grid_vector().size();

            typename direct_luts::C_aoa_long_aos_dh_mid::table_type table(size_x * size_y * size_z);
            auto t = std::begin(table);
            for (auto dh : { "dhm25", "dhm10", "dh0", "dh10", "dh25" }) {
                const auto [CX, num_rows_CX, num_columns_CX] = read_csv_table(path_dataset + "/CX_aoa_aos_" + dh + ".csv");
                const auto [CZ, num_rows_CZ, num_columns_CZ] = read_csv_table(path_dataset + "/CZ_aoa_aos_" + dh + ".csv");
                const auto [CM, num_rows_CM, num_columns_CM] = read_csv_table(path_dataset + "/CM_aoa_aos_" + dh + ".csv");

                if (num_rows_CX != size_x || num_rows_CZ != size_x || num_rows_CM != size_x ||
                    num_columns_CX != size_y || num_columns_CZ != size_y || num_columns_CM != size_y) {

                    TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                            " reading \"direct_luts::C_aoa_long_aos_dh_mid\".")

                }

                for (auto count = 0u; count < size_x * size_y; ++count) {
                    *t++ = { CX[count], CZ[count], CM[count] };
                }

            }

            _dluts.C_aoa_long_aos_dh_mid.build(std::move(table), size_x, size_y, size_z);

        }


        // read tables for "C_aoa_long_aos_dh_short"
        {
            const auto size_x = _prelookups.aoa_long_deg.grid_vector().size();
            const auto size_y = _prelookups.aos_deg.grid_vector().size();
            const auto size_z = _prelookups.dh_short_deg.grid_vector().size();

            typename direct_luts::C_aoa_long_aos_dh_short::table_type table(size_x * size_y * size_z);
            auto t = std::begin(table);
            for (auto dh : { "dhm25", "dh0", "dh25" }) {
                const auto [CLL, num_rows_CLL, num_columns_CLL] = read_csv_table(path_dataset + "/CLL_aoa_aos_" + dh + ".csv");
                const auto [CN, num_rows_CN, num_columns_CN] = read_csv_table(path_dataset + "/CN_aoa_aos_" + dh + ".csv");

                if (num_rows_CLL != size_x || num_rows_CLL != size_x ||
                    num_columns_CN != size_y || num_columns_CN != size_y) {

                    TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                            " reading \"direct_luts::C_aoa_long_aos_dh_short\".")

                }

                for (auto count = 0u; count < size_x * size_y; ++count) {
                    *t++ = { CLL[count], CN[count] };
                }

            }

            _dluts.C_aoa_long_aos_dh_short.build(std::move(table), size_x, size_y, size_z);

        }


        // read tables for "C_aoa_long_aos"
        {
            const auto size_x = _prelookups.aoa_long_deg.grid_vector().size();
            const auto size_y = _prelookups.aos_deg.grid_vector().size();

            const auto [CXdh0, num_rows_CXdh0, num_columns_CXdh0] = read_csv_table(path_dataset + "/CX_aoa_aos_dh0.csv");
            const auto [CY, num_rows_CY, num_columns_CY] = read_csv_table(path_dataset + "/CY_aoa_aos.csv");
            const auto [CYda20, num_rows_CYda20, num_columns_CYda20] = read_csv_table(path_dataset + "/CYda20_aoa_aos.csv");
            const auto [CYdr30, num_rows_CYdr30, num_columns_CYdr30] = read_csv_table(path_dataset + "/CYdr30_aoa_aos.csv");
            const auto [CZdh0, num_rows_CZdh0, num_columns_CZdh0] = read_csv_table(path_dataset + "/CZ_aoa_aos_dh0.csv");
            const auto [CLLdh0, num_rows_CLLdh0, num_columns_CLLdh0] = read_csv_table(path_dataset + "/CLL_aoa_aos_dh0.csv");
            const auto [CLLda20, num_rows_CLLda20, num_columns_CLLda20] = read_csv_table(path_dataset + "/CLLda20_aoa_aos.csv");
            const auto [CLLdr30, num_rows_CLLdr30, num_columns_CLLdr30] = read_csv_table(path_dataset + "/CLLdr30_aoa_aos.csv");
            const auto [CMdh0, num_rows_CMdh0, num_columns_CMdh0] = read_csv_table(path_dataset + "/CM_aoa_aos_dh0.csv");
            const auto [CNdh0, num_rows_CNdh0, num_columns_CNdh0] = read_csv_table(path_dataset + "/CN_aoa_aos_dh0.csv");
            const auto [CNda20, num_rows_CNda20, num_columns_CNda20] = read_csv_table(path_dataset + "/CNda20_aoa_aos.csv");
            const auto [CNdr30, num_rows_CNdr30, num_columns_CNdr30] = read_csv_table(path_dataset + "/CNdr30_aoa_aos.csv");

            if (num_rows_CXdh0 != size_x ||
                num_rows_CY != size_x || num_rows_CYda20 != size_x || num_rows_CYdr30 != size_x ||
                num_rows_CZdh0 != size_x ||
                num_rows_CLLdh0 != size_x || num_rows_CLLda20 != size_x || num_rows_CLLdr30 != size_x ||
                num_rows_CMdh0 != size_x ||
                num_rows_CNdh0 != size_x || num_rows_CNda20 != size_x || num_rows_CNdr30 != size_x ||
                num_columns_CXdh0 != size_y ||
                num_columns_CY != size_y || num_columns_CYda20 != size_y || num_columns_CYdr30 != size_y ||
                num_columns_CZdh0 != size_y ||
                num_columns_CLLdh0 != size_y || num_columns_CLLda20 != size_y || num_columns_CLLdr30 != size_y ||
                num_columns_CMdh0 != size_y ||
                num_columns_CNdh0 != size_y || num_columns_CNda20 != size_y || num_columns_CNdr30 != size_y) {

                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                        " reading \"direct_luts::C_aoa_long_aos\".")

            }

            typename direct_luts::C_aoa_long_aos::table_type t(size_x * size_y);
            for (auto count = 0u; count < size_x * size_y; ++count) {
                t[count] = { CXdh0[count],
                             CY[count], CYda20[count], CYdr30[count],
                             CZdh0[count],
                             CLLdh0[count], CLLda20[count], CLLdr30[count],
                             CMdh0[count],
                             CNdh0[count], CNda20[count], CNdr30[count] };

            }

            _dluts.C_aoa_long_aos.build(std::move(t), size_x, size_y);

        }


        // read tables for "C_aoa_short_aos"
        {
            const auto size_x = _prelookups.aoa_short_deg.grid_vector().size();
            const auto size_y = _prelookups.aos_deg.grid_vector().size();

            const auto [CXlef, num_rows_CXlef, num_columns_CXlef] = read_csv_table(path_dataset + "/CXlef_aoa_aos.csv");
            const auto [CYlef, num_rows_CYlef, num_columns_CYlef] = read_csv_table(path_dataset + "/CYlef_aoa_aos.csv");
            const auto [CYda20lef, num_rows_CYda20lef, num_columns_CYda20lef] = read_csv_table(path_dataset + "/CYda20lef_aoa_aos.csv");
            const auto [CZlef, num_rows_CZlef, num_columns_CZlef] = read_csv_table(path_dataset + "/CZlef_aoa_aos.csv");
            const auto [CLLlef, num_rows_CLLlef, num_columns_CLLlef] = read_csv_table(path_dataset + "/CLLlef_aoa_aos.csv");
            const auto [CLLda20lef, num_rows_CLLda20lef, num_columns_CLLda20lef] = read_csv_table(path_dataset + "/CLLda20lef_aoa_aos.csv");
            const auto [CMlef, num_rows_CMlef, num_columns_CMlef] = read_csv_table(path_dataset + "/CMlef_aoa_aos.csv");
            const auto [CNlef, num_rows_CNlef, num_columns_CNlef] = read_csv_table(path_dataset + "/CNlef_aoa_aos.csv");
            const auto [CNda20lef, num_rows_CNda20lef, num_columns_CNda20lef] = read_csv_table(path_dataset + "/CNda20lef_aoa_aos.csv");

            if (num_rows_CXlef != size_x ||
                num_rows_CYlef != size_x || num_rows_CYda20lef != size_x ||
                num_rows_CZlef != size_x ||
                num_rows_CLLlef != size_x || num_rows_CLLda20lef != size_x ||
                num_rows_CMlef != size_x ||
                num_rows_CNlef != size_x || num_rows_CNda20lef != size_x ||
                num_columns_CXlef != size_y ||
                num_columns_CYlef != size_y || num_columns_CYda20lef != size_y ||
                num_columns_CZlef != size_y ||
                num_columns_CLLlef != size_y || num_columns_CLLda20lef != size_y ||
                num_columns_CMlef != size_y ||
                num_columns_CNlef != size_y || num_columns_CNda20lef != size_y) {

                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                        " reading \"direct_luts::C_aoa_short_aos\".")

            }

            typename direct_luts::C_aoa_short_aos::table_type t(size_x * size_y);
            for (auto count = 0u; count < size_x * size_y; ++count) {
                t[count] = { CXlef[count],
                             CYlef[count], CYda20lef[count],
                             CZlef[count],
                             CLLlef[count], CLLda20lef[count],
                             CMlef[count],
                             CNlef[count], CNda20lef[count] };

            }

            _dluts.C_aoa_short_aos.build(std::move(t), size_x, size_y);

        }


        // read tables for "C_aoa_long_dh_long"
        {
            const auto size_x = _prelookups.aoa_long_deg.grid_vector().size();
            const auto size_y = _prelookups.dh_long_deg.grid_vector().size();

            const auto [DCMds, num_rows_DCMds, num_columns_DCMds] = read_csv_table(path_dataset + "/DCMds_aoa_dh.csv");

            if (num_rows_DCMds != size_x || num_columns_DCMds != size_y) {
                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                        " reading \"direct_luts::C_aoa_long_dh_long\".")

            }

            typename direct_luts::C_aoa_long_dh_long::table_type t(size_x * size_y);
            std::transform(std::cbegin(DCMds), std::cend(DCMds), std::begin(t),
                          [](auto val) { return typename direct_luts::C_aoa_long_dh_long::table_type::value_type{ val }; });

            _dluts.C_aoa_long_dh_long.build(std::move(t), size_x, size_y);

        }


        // read tables for "C_aoa_long"
        {
            const auto size_x = _prelookups.aoa_long_deg.grid_vector().size();
            constexpr auto size_y = std::size_t{ 1u };

            const auto [CXq, num_rows_CXq, num_columns_CXq] = read_csv_table(path_dataset + "/CXq_aoa.csv");
            const auto [DCXsb, num_rows_DCXsb, num_columns_DCXsb] = read_csv_table(path_dataset + "/DCXsb_aoa.csv");
            const auto [CYp, num_rows_CYp, num_columns_CYp] = read_csv_table(path_dataset + "/CYp_aoa.csv");
            const auto [CYr, num_rows_CYr, num_columns_CYr] = read_csv_table(path_dataset + "/CYr_aoa.csv");
            const auto [CZq, num_rows_CZq, num_columns_CZq] = read_csv_table(path_dataset + "/CZq_aoa.csv");
            const auto [DCZsb, num_rows_DCZsb, num_columns_DCZsb] = read_csv_table(path_dataset + "/DCZsb_aoa.csv");
            const auto [CLLp, num_rows_CLLp, num_columns_CLLp] = read_csv_table(path_dataset + "/CLLp_aoa.csv");
            const auto [CLLr, num_rows_CLLr, num_columns_CLLr] = read_csv_table(path_dataset + "/CLLr_aoa.csv");
            const auto [DCLLbeta, num_rows_DCLLbeta, num_columns_DCLLbeta] = read_csv_table(path_dataset + "/DCLLbeta_aoa.csv");
            const auto [CMq, num_rows_CMq, num_columns_CMq] = read_csv_table(path_dataset + "/CMq_aoa.csv");
            const auto [DCM, num_rows_DCM, num_columns_DCM] = read_csv_table(path_dataset + "/DCM_aoa.csv");
            const auto [DCMsb, num_rows_DCMsb, num_columns_DCMsb] = read_csv_table(path_dataset + "/DCMsb_aoa.csv");
            const auto [CNp, num_rows_CNp, num_columns_CNp] = read_csv_table(path_dataset + "/CNp_aoa.csv");
            const auto [CNr, num_rows_CNr, num_columns_CNr] = read_csv_table(path_dataset + "/CNr_aoa.csv");
            const auto [DCNbeta, num_rows_DCNbeta, num_columns_DCNbeta] = read_csv_table(path_dataset + "/DCNbeta_aoa.csv");
            const auto [DCNda, num_rows_DCNda, num_columns_DCNda] = read_csv_table(path_dataset + "/DCNda_aoa.csv");

            if (num_rows_CXq != size_x || num_rows_DCXsb != size_x ||
                num_rows_CYp != size_x || num_rows_CYr != size_x ||
                num_rows_CZq != size_x || num_rows_DCZsb != size_x ||
                num_rows_CLLp != size_x || num_rows_CLLr != size_x || num_rows_DCLLbeta != size_x ||
                num_rows_CMq != size_x || num_rows_DCM != size_x || num_rows_DCMsb != size_x ||
                num_rows_CNp != size_x || num_rows_CNr != size_x || num_rows_DCNbeta != size_x || num_rows_DCNda != size_x ||
                num_columns_CXq != size_y || num_columns_DCXsb != size_y ||
                num_columns_CYp != size_y || num_columns_CYr != size_y ||
                num_columns_CZq != size_y || num_columns_DCZsb != size_y ||
                num_columns_CLLp != size_y || num_columns_CLLr != size_y || num_columns_DCLLbeta != size_y ||
                num_columns_CMq != size_y || num_columns_DCM != size_y || num_columns_DCMsb != size_y ||
                num_columns_CNp != size_y || num_columns_CNr != size_y || num_columns_DCNbeta != size_y || num_columns_DCNda != size_y) {

                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                        " reading \"direct_luts::C_aoa_long\".")

            }

            typename direct_luts::C_aoa_long::table_type t(size_x);
            for (auto count = 0u; count < size_x; ++count) {
                t[count] = { CXq[count], DCXsb[count],
                             CYp[count], CYr[count],
                             CZq[count], DCZsb[count],
                             CLLp[count], CLLr[count], DCLLbeta[count],
                             CMq[count], DCM[count], DCMsb[count],
                             CNp[count], CNr[count], DCNbeta[count], DCNda[count] };

            }

            _dluts.C_aoa_long.build(std::move(t), size_x);

        }


        // read tables for "C_aoa_short"
        {
            const auto size_x = _prelookups.aoa_short_deg.grid_vector().size();
            constexpr auto size_y = std::size_t{ 1u };

            const auto [DCXqlef, num_rows_DCXqlef, num_columns_DCXqlef] = read_csv_table(path_dataset + "/DCXqlef_aoa.csv");
            const auto [DCYplef, num_rows_DCYplef, num_columns_DCYplef] = read_csv_table(path_dataset + "/DCYplef_aoa.csv");
            const auto [DCYrlef, num_rows_DCYrlef, num_columns_DCYrlef] = read_csv_table(path_dataset + "/DCYrlef_aoa.csv");
            const auto [DCZqlef, num_rows_DCZqlef, num_columns_DCZqlef] = read_csv_table(path_dataset + "/DCZqlef_aoa.csv");
            const auto [DCLLplef, num_rows_DCLLplef, num_columns_DCLLplef] = read_csv_table(path_dataset + "/DCLLplef_aoa.csv");
            const auto [DCLLrlef, num_rows_DCLLrlef, num_columns_DCLLrlef] = read_csv_table(path_dataset + "/DCLLrlef_aoa.csv");
            const auto [DCMqlef, num_rows_DCMqlef, num_columns_DCMqlef] = read_csv_table(path_dataset + "/DCMqlef_aoa.csv");
            const auto [DCNplef, num_rows_DCNplef, num_columns_DCNplef] = read_csv_table(path_dataset + "/DCNplef_aoa.csv");
            const auto [DCNrlef, num_rows_DCNrlef, num_columns_DCNrlef] = read_csv_table(path_dataset + "/DCNrlef_aoa.csv");
            const auto [CLalphadot, num_rows_CLalphadot, num_columns_CLalphadot] = read_csv_table(path_dataset + "/CLalphadot_aoa.csv");
            const auto [CMalphadot, num_rows_CMalphadot, num_columns_CMalphadot] = read_csv_table(path_dataset + "/CMalphadot_aoa.csv");

            if (num_rows_DCXqlef != size_x ||
                num_rows_DCYplef != size_x || num_rows_DCYrlef != size_x ||
                num_rows_DCZqlef != size_x ||
                num_rows_DCLLplef != size_x || num_rows_DCLLrlef != size_x ||
                num_rows_DCMqlef != size_x ||
                num_rows_DCNplef != size_x || num_rows_DCNrlef != size_x ||
                num_rows_CLalphadot != size_x || num_rows_CMalphadot != size_x ||
                num_columns_DCXqlef != size_y ||
                num_columns_DCYplef != size_y || num_columns_DCYrlef != size_y ||
                num_columns_DCZqlef != size_y ||
                num_columns_DCLLplef != size_y || num_columns_DCLLrlef != size_y ||
                num_columns_DCMqlef != size_y ||
                num_columns_DCNplef != size_y || num_columns_DCNrlef != size_y ||
                num_columns_CLalphadot != size_y || num_columns_CMalphadot != size_y) {

                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                        " reading \"direct_luts::C_aoa_short\".")

            }

            typename direct_luts::C_aoa_short::table_type t(size_x);
            for (auto count = 0u; count < size_x; ++count) {
                t[count] = { DCXqlef[count],
                             DCYplef[count], DCYrlef[count],
                             DCZqlef[count],
                             DCLLplef[count], DCLLrlef[count],
                             DCMqlef[count],
                             DCNplef[count], DCNrlef[count],
                             CLalphadot[count], CMalphadot[count] };

            }

            _dluts.C_aoa_short.build(std::move(t), size_x);

        }


        // read tables for "C_dh_mid"
        {
            const auto size_x = _prelookups.dh_mid_deg.grid_vector().size();
            constexpr auto size_y  = std::size_t{ 1u };

            const auto [etadh, num_rows_etadh, num_columns_etadh] = read_csv_table(path_dataset + "/etadh_dh.csv");

            if (num_rows_etadh != size_x || num_columns_etadh != size_y) {
                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_aerodynamics::build: inconsistent sizes while"
                                        " reading \"direct_luts::C_dh_mid\".")

            }

            typename direct_luts::C_dh_mid::table_type t(size_x);
            std::transform(std::cbegin(etadh), std::cend(etadh), std::begin(t),
                           [](auto val) { return typename direct_luts::C_dh_mid::table_type::value_type{ val }; });

            _dluts.C_dh_mid.build(std::move(t), size_x);

        }

    }


    template<typename T>
    constexpr auto F16_Nguyen_aerodynamics::outputs(const T &dynamic_pressure_Pa,
                                                    const T &TAS_mps, const T &aoa_deg, const T &aos_deg,
                                                    const T &paero_radps, const T &qaero_radps, const T &raero_radps,
                                                    const T &dh_deg, const T &dlef_deg, const T &dsb_deg,
                                                    const T &da_deg, const T &dr_deg,
                                                    const T &aoadot_radps,
                                                    const T &xcg_per_MAC,
                                                    const T &ycg_per_semiwingspan,
                                                    const T &MAC_m,
                                                    const T &wingspan_m,
                                                    const T &wing_surface_m2) const -> outputs_type<T>
    {
        //
        // Evaluates the model's output equation, i.e.,
        // "y = G(x, u, ...)", and returns the model's "y".
        //

        // evaluate the total force & moment coefficients
        const auto [CXtot, CYtot, CZtot,
                    CLLtot, CMtot, CNtot] = force_and_moment_coefficients(aoa_deg, aos_deg,
                                                                          dh_deg,
                                                                          TAS_mps,
                                                                          paero_radps, qaero_radps, raero_radps,
                                                                          dlef_deg, dsb_deg,
                                                                          da_deg, dr_deg,
                                                                          aoadot_radps,
                                                                          xcg_per_MAC,
                                                                          ycg_per_semiwingspan,
                                                                          MAC_m,
                                                                          wingspan_m);


        // calculate the forces & moments, return the model's outputs
        const auto qS = dynamic_pressure_Pa * wing_surface_m2;

        return outputs_type<T>{ CXtot, CYtot, CZtot,
                                CLLtot, CMtot, CNtot,
                                qS * CXtot, qS * CYtot, qS * CZtot,
                                qS * wingspan_m * CLLtot, qS * MAC_m * CMtot, qS * wingspan_m * CNtot };

    }


    template<typename T>
    constexpr auto F16_Nguyen_aerodynamics::force_and_moment_coefficients(const T &aoa_deg, const T &aos_deg,
                                                                          const T &dh_deg,
                                                                          const T &TAS_mps,
                                                                          const T &paero_radps, const T &qaero_radps, const T &raero_radps,
                                                                          const T &dlef_deg, const T &dsb_deg,
                                                                          const T &da_deg, const T &dr_deg,
                                                                          const T &aoadot_radps,
                                                                          const T &xcg_per_MAC,
                                                                          const T &ycg_per_semiwingspan,
                                                                          const T &MAC_m,
                                                                          const T &wingspan_m) const -> force_and_moment_coefficients_type<T>
    {
        //
        // Evaluates the model's total force & moment coefficients,
        // and returns them as an std::array.
        //

        // calculate the nondimensional surfaces, xcg and rates
        const auto dlef_nondim = T{ 1 } - dlef_deg / T{ default_parameters::reference_leading_edge_flap_deg };
        const auto dsb_nondim = dsb_deg / T{ default_parameters::reference_speedbreak_deg };
        const auto da_nondim = da_deg / T{ default_parameters::reference_aileron_deg };
        const auto dr_nondim = dr_deg / T{ default_parameters::reference_rudder_deg };

        const auto Dxcg_nondim = T{ default_parameters::reference_xcg_per_MAC } - xcg_per_MAC;
        const auto Dycg_nondim = T{ 0.5 } * (T{ default_parameters::reference_ycg_per_semiwingspan } - ycg_per_semiwingspan);

        const auto mid_invTAS = T{ 0.5 } / TAS_mps;
        const auto p_nondim = paero_radps * wingspan_m * mid_invTAS;
        const auto q_nondim = qaero_radps * MAC_m * mid_invTAS;
        const auto r_nondim = raero_radps * wingspan_m * mid_invTAS;
        const auto aoadot_nondim = aoadot_radps * MAC_m * mid_invTAS;


        // evaluate all the dataset's coefficients
        const auto c = dataset_coefficients(aoa_deg, aos_deg,
                                            dh_deg);


        // calculate the total force and moment coefficients
        const auto CXtot = c[dataset_coefficient_names::CX] +
                           (c[dataset_coefficient_names::CXlef] - c[dataset_coefficient_names::CXdh0]) * dlef_nondim +
                           c[dataset_coefficient_names::DCXsb] * dsb_nondim +
                           (c[dataset_coefficient_names::CXq] + c[dataset_coefficient_names::DCXqlef] * dlef_nondim) * q_nondim +
                           c[dataset_coefficient_names::CLalphadot] * tr7::AD_sin(tr7::deg2rad(aoa_deg)) * aoadot_nondim;

        const auto CZtot = c[dataset_coefficient_names::CZ] +
                           (c[dataset_coefficient_names::CZlef] - c[dataset_coefficient_names::CZdh0]) * dlef_nondim +
                           c[dataset_coefficient_names::DCZsb] * dsb_nondim +
                           (c[dataset_coefficient_names::CZq] + c[dataset_coefficient_names::DCZqlef] * dlef_nondim) * q_nondim -
                           c[dataset_coefficient_names::CLalphadot] * tr7::AD_cos(tr7::deg2rad(aoa_deg)) * aoadot_nondim;

        const auto CMtot = c[dataset_coefficient_names::CM] * c[dataset_coefficient_names::etadh] +
                           CZtot * Dxcg_nondim +
                           (c[dataset_coefficient_names::CMlef] - c[dataset_coefficient_names::CMdh0]) * dlef_nondim +
                           c[dataset_coefficient_names::DCMsb] * dsb_nondim +
                           (c[dataset_coefficient_names::CMq] + c[dataset_coefficient_names::DCMqlef] * dlef_nondim) * q_nondim +
                           c[dataset_coefficient_names::DCM] +
                           c[dataset_coefficient_names::DCMds] +
                           c[dataset_coefficient_names::CMalphadot] * aoadot_nondim;

        const auto DCYda20 = c[dataset_coefficient_names::CYda20] - c[dataset_coefficient_names::CY];
        const auto CYtot = c[dataset_coefficient_names::CY] +
                           (c[dataset_coefficient_names::CYlef] - c[dataset_coefficient_names::CY]) * dlef_nondim +
                           (DCYda20 + (c[dataset_coefficient_names::CYda20lef] - c[dataset_coefficient_names::CYlef] - DCYda20) * dlef_nondim) * da_nondim +
                           (c[dataset_coefficient_names::CYdr30] - c[dataset_coefficient_names::CY]) * dr_nondim +
                           (c[dataset_coefficient_names::CYr] + c[dataset_coefficient_names::DCYrlef] * dlef_nondim) * r_nondim +
                           (c[dataset_coefficient_names::CYp] + c[dataset_coefficient_names::DCYplef] * dlef_nondim) * p_nondim;

        const auto DCNda20 = c[dataset_coefficient_names::CNda20] - c[dataset_coefficient_names::CNdh0];
        const auto CNtot = c[dataset_coefficient_names::CN] +
                           CXtot * Dycg_nondim - CYtot * Dxcg_nondim * MAC_m / wingspan_m +
                           (c[dataset_coefficient_names::CNlef] - c[dataset_coefficient_names::CNdh0]) * dlef_nondim +
                           (DCNda20 + (c[dataset_coefficient_names::CNda20lef] - c[dataset_coefficient_names::CNlef] - DCNda20) * dlef_nondim) * da_nondim +
                           (c[dataset_coefficient_names::CNdr30] - c[dataset_coefficient_names::CNdh0]) * dr_nondim +
                           (c[dataset_coefficient_names::CNr] + c[dataset_coefficient_names::DCNrlef] * dlef_nondim) * r_nondim +
                           (c[dataset_coefficient_names::CNp] + c[dataset_coefficient_names::DCNplef] * dlef_nondim) * p_nondim +
                           c[dataset_coefficient_names::DCNbeta] * aos_deg;

        const auto DCLLda20 = c[dataset_coefficient_names::CLLda20] - c[dataset_coefficient_names::CLLdh0];
        const auto CLLtot = c[dataset_coefficient_names::CLL] -
                            CZtot * Dycg_nondim +
                            (c[dataset_coefficient_names::CLLlef] - c[dataset_coefficient_names::CLLdh0]) * dlef_nondim +
                            (DCLLda20 + (c[dataset_coefficient_names::CLLda20lef] - c[dataset_coefficient_names::CLLlef] - DCLLda20) * dlef_nondim) * da_nondim +
                            (c[dataset_coefficient_names::CLLdr30] - c[dataset_coefficient_names::CLLdh0]) * dr_nondim +
                            (c[dataset_coefficient_names::CLLr] + c[dataset_coefficient_names::DCLLrlef] * dlef_nondim) * r_nondim +
                            (c[dataset_coefficient_names::CLLp] + c[dataset_coefficient_names::DCLLplef] * dlef_nondim) * p_nondim +
                            c[dataset_coefficient_names::DCLLbeta] * aos_deg;


        // return them
        return force_and_moment_coefficients_type<T>{ CXtot, CYtot, CZtot,
                                                      CLLtot, CMtot, CNtot };

    }


    template<typename T>
    constexpr auto F16_Nguyen_aerodynamics::dataset_coefficients(const T &aoa_deg, const T &aos_deg,
                                                                 const T &dh_deg) const -> dataset_coefficients_type<T>
    {
        //
        // Evaluates the dataset's coefficients, returning a
        // central difference approximation when we're in AD.
        //

        if constexpr (!tr7::is_AD_scalar<T>::value) {
            return lookup_dataset_coefficients(aoa_deg, aos_deg, dh_deg);

        }
        else {
            // the AD with a linear lookup table behaves as a forward difference,
            // we need to add the influence of the backward one to get a central
            // result
            auto c = lookup_dataset_coefficients(aoa_deg, aos_deg, dh_deg);

            const auto aoa_numerical_deg = tr7::AD_force_scalar_value(aoa_deg);
            const auto aos_numerical_deg = tr7::AD_force_scalar_value(aos_deg);
            const auto dh_numerical_deg = tr7::AD_force_scalar_value(dh_deg);

            constexpr auto delta_numjac = decltype(aoa_numerical_deg){ 1e-9 };
            const auto c_aoa_m = lookup_dataset_coefficients(aoa_numerical_deg - delta_numjac, aos_numerical_deg, dh_numerical_deg);
            const auto c_aos_m = lookup_dataset_coefficients(aoa_numerical_deg, aos_numerical_deg - delta_numjac, dh_numerical_deg);
            const auto c_dh_m = lookup_dataset_coefficients(aoa_numerical_deg, aos_numerical_deg, dh_numerical_deg - delta_numjac);

            const auto Daoa_div_delta = (aoa_deg - aoa_numerical_deg) / delta_numjac;
            const auto Daos_div_delta = (aos_deg - aos_numerical_deg) / delta_numjac;
            const auto Ddh_div_delta = (dh_deg - dh_numerical_deg) / delta_numjac;

            tr7::tuple_index_for(c, [&](auto &ci, auto i)
                                       {
                                           const auto ci_numerical = tr7::AD_force_scalar_value(ci);
                                           const auto Dc_aoa = ci_numerical - c_aoa_m[i];
                                           const auto Dc_aos = ci_numerical - c_aos_m[i];
                                           const auto Dc_dh = ci_numerical - c_dh_m[i];

                                           ci = T{ 0.5 } * ci +
                                                decltype(ci_numerical){ 0.5 } * (ci_numerical +
                                                                                 Dc_aoa * Daoa_div_delta +
                                                                                 Dc_aos * Daos_div_delta +
                                                                                 Dc_dh * Ddh_div_delta);

                                       });

            return c;

        }

    }


    template<typename T>
    constexpr auto F16_Nguyen_aerodynamics::lookup_dataset_coefficients(const T &aoa_deg, const T &aos_deg,
                                                                        const T &dh_deg) const -> dataset_coefficients_type<T>
    {
        //
        // Evaluates the dataset's coefficients and returns
        // an std::array holding all of them.
        //

        // query prelookups
        const auto p_aoa_long = _prelookups.aoa_long_deg.template evaluate<true>(aoa_deg);
        const auto p_aoa_short = std::get<1>(p_aoa_long.first) != tr7::lookup::extrapolation_flag::upper ? p_aoa_long :
                                                                                                           _prelookups.aoa_short_deg.template evaluate<true>(aoa_deg);
        const auto p_aos = _prelookups.aos_deg.template evaluate<true>(aos_deg);
        const auto p_dh_long = _prelookups.dh_long_deg.template evaluate<true>(dh_deg);
        const auto p_dh_mid = _prelookups.dh_mid_deg.template evaluate<true>(dh_deg);
        const auto p_dh_short = _prelookups.dh_short_deg.template evaluate<true>(dh_deg);


        // query luts
        const auto C_aoa_long_aos_dh_mid = _dluts.C_aoa_long_aos_dh_mid(p_aoa_long, p_aos, p_dh_mid);
        const auto C_aoa_long_aos_dh_short = _dluts.C_aoa_long_aos_dh_short(p_aoa_long, p_aos, p_dh_short);
        const auto C_aoa_long_aos = _dluts.C_aoa_long_aos(p_aoa_long, p_aos);
        const auto C_aoa_short_aos = _dluts.C_aoa_short_aos(p_aoa_short, p_aos);
        const auto C_aoa_long_dh_long = _dluts.C_aoa_long_dh_long(p_aoa_long, p_dh_long);
        const auto C_aoa_long = _dluts.C_aoa_long(p_aoa_long);
        const auto C_aoa_short = _dluts.C_aoa_short(p_aoa_short);
        const auto C_dh_mid = _dluts.C_dh_mid(p_dh_mid);


        // get all coefficients and return them
#define GET_CAST_COEFF(C, NAME) T{ C[direct_luts::C::NAME] }
        return dataset_coefficients_type<T>{ // CX-related
                                             GET_CAST_COEFF(C_aoa_long_aos_dh_mid, CX),
                                             GET_CAST_COEFF(C_aoa_long_aos, CXdh0),
                                             GET_CAST_COEFF(C_aoa_short_aos, CXlef),
                                             GET_CAST_COEFF(C_aoa_long, CXq),
                                             GET_CAST_COEFF(C_aoa_long, DCXsb),
                                             GET_CAST_COEFF(C_aoa_short, DCXqlef),

                                             // CY-related
                                             GET_CAST_COEFF(C_aoa_long_aos, CY),
                                             GET_CAST_COEFF(C_aoa_long_aos, CYda20),
                                             GET_CAST_COEFF(C_aoa_long_aos, CYdr30),
                                             GET_CAST_COEFF(C_aoa_short_aos, CYlef),
                                             GET_CAST_COEFF(C_aoa_short_aos, CYda20lef),
                                             GET_CAST_COEFF(C_aoa_long, CYp),
                                             GET_CAST_COEFF(C_aoa_long, CYr),
                                             GET_CAST_COEFF(C_aoa_short, DCYplef),
                                             GET_CAST_COEFF(C_aoa_short, DCYrlef),

                                             // CZ-related
                                             GET_CAST_COEFF(C_aoa_long_aos_dh_mid, CZ),
                                             GET_CAST_COEFF(C_aoa_long_aos, CZdh0),
                                             GET_CAST_COEFF(C_aoa_short_aos, CZlef),
                                             GET_CAST_COEFF(C_aoa_long, CZq),
                                             GET_CAST_COEFF(C_aoa_long, DCZsb),
                                             GET_CAST_COEFF(C_aoa_short, DCZqlef),

                                             // CLL-related
                                             GET_CAST_COEFF(C_aoa_long_aos_dh_short, CLL),
                                             GET_CAST_COEFF(C_aoa_long_aos, CLLdh0),
                                             GET_CAST_COEFF(C_aoa_long_aos, CLLda20),
                                             GET_CAST_COEFF(C_aoa_long_aos, CLLdr30),
                                             GET_CAST_COEFF(C_aoa_short_aos, CLLlef),
                                             GET_CAST_COEFF(C_aoa_short_aos, CLLda20lef),
                                             GET_CAST_COEFF(C_aoa_long, CLLp),
                                             GET_CAST_COEFF(C_aoa_long, CLLr),
                                             GET_CAST_COEFF(C_aoa_long, DCLLbeta),
                                             GET_CAST_COEFF(C_aoa_short, DCLLplef),
                                             GET_CAST_COEFF(C_aoa_short, DCLLrlef),

                                             // CM-related
                                             GET_CAST_COEFF(C_aoa_long_aos_dh_mid, CM),
                                             GET_CAST_COEFF(C_aoa_long_aos, CMdh0),
                                             GET_CAST_COEFF(C_aoa_short_aos, CMlef),
                                             GET_CAST_COEFF(C_aoa_long_dh_long, DCMds),
                                             GET_CAST_COEFF(C_aoa_long, CMq),
                                             GET_CAST_COEFF(C_aoa_long, DCM),
                                             GET_CAST_COEFF(C_aoa_long, DCMsb),
                                             GET_CAST_COEFF(C_aoa_short, DCMqlef),

                                             // CN-related
                                             GET_CAST_COEFF(C_aoa_long_aos_dh_short, CN),
                                             GET_CAST_COEFF(C_aoa_long_aos, CNdh0),
                                             GET_CAST_COEFF(C_aoa_long_aos, CNda20),
                                             GET_CAST_COEFF(C_aoa_long_aos, CNdr30),
                                             GET_CAST_COEFF(C_aoa_short_aos, CNlef),
                                             GET_CAST_COEFF(C_aoa_short_aos, CNda20lef),
                                             GET_CAST_COEFF(C_aoa_long, CNp),
                                             GET_CAST_COEFF(C_aoa_long, CNr),
                                             GET_CAST_COEFF(C_aoa_long, DCNbeta),
                                             GET_CAST_COEFF(C_aoa_long, DCNda),
                                             GET_CAST_COEFF(C_aoa_short, DCNplef),
                                             GET_CAST_COEFF(C_aoa_short, DCNrlef),

                                             // etadh
                                             GET_CAST_COEFF(C_dh_mid, etadh),

                                             // aoadot ones
                                             GET_CAST_COEFF(C_aoa_short, CLalphadot),
                                             GET_CAST_COEFF(C_aoa_short, CMalphadot) };
    // macro cleanup
#undef GET_CAST_COEFF

    }


}


#endif