#include <iomanip>

#include "F16_Nguyen/F16_Nguyen_plant.h"


/*
int main()
{
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);

    const std::string aeropath = "/home/c83833/codes/F16/F16_Nguyen/datasets/aero";
    F16_Nguyen::F16_Nguyen_aerodynamics aero(aeropath);

    // validate CX CM CZ
    {
        const auto alphas = std::get<0>(F16_Nguyen::read_csv_table(aeropath + "/aoa_long_breaks_deg.csv"));
        const auto betas = std::get<0>(F16_Nguyen::read_csv_table(aeropath + "/aos_breaks_deg.csv"));
        const auto dhs = std::get<0>(F16_Nguyen::read_csv_table(aeropath + "/dh_mid_breaks_deg.csv"));
        const auto num_rows = alphas.size();
        const auto num_cols = betas.size();
        auto dh = std::begin(dhs);
        for (auto dh_str : { "dhm25", "dhm10", "dh0", "dh10", "dh25" }) {
            const auto CX = std::get<0>(F16_Nguyen::read_csv_table(aeropath + "/CX_aoa_aos_" + dh_str + ".csv"));
            const auto CZ = std::get<0>(F16_Nguyen::read_csv_table(aeropath + "/CZ_aoa_aos_" + dh_str + ".csv"));
            const auto CM = std::get<0>(F16_Nguyen::read_csv_table(aeropath + "/CM_aoa_aos_" + dh_str + ".csv"));
        
            for (auto col = 0u; col != num_cols; ++col) {
                for (auto row = 0u; row != num_rows; ++row) {
                    const auto coeffs = aero.dataset(alphas[row], betas[col], *dh);
                    const auto ij = row + num_rows * col;
        
                    if (!tr7::eq_fp(coeffs[F16_Nguyen::F16_Nguyen_aerodynamics<>::dataset_coefficients::CX], CX[ij], 10.)) {
                        TR7_THROW_RUNTIME_ERROR("MAL CX");
                    }
                    if (!tr7::eq_fp(coeffs[F16_Nguyen::F16_Nguyen_aerodynamics<>::dataset_coefficients::CZ], CZ[ij], 10.)) {
                        TR7_THROW_RUNTIME_ERROR("MAL CZ");
                    }
                    if (!tr7::eq_fp(coeffs[F16_Nguyen::F16_Nguyen_aerodynamics<>::dataset_coefficients::CM], CM[ij], 10.)) {
                        TR7_THROW_RUNTIME_ERROR("MAL CM");
                    }
        
                }
            }
        
            ++dh;
        
        }

    }


    std::cout << "done." << '\n';


}
*/


int main()
{
    const std::string aeropath = "../../datasets/aero";
    const std::string enginepath = "../../datasets/engine";

    F16_Nguyen::F16_Nguyen_plant pl(aeropath, enginepath);


    auto sl = pl.state_names2str().cbegin();
    std::cout << "plant_num_states  = " << decltype(pl)::num_states() << ";\n";
    std::cout << "plant_state_names = {'" << *sl++ << "'";
    while (sl != pl.state_names2str().cend()) {
        std::cout << "; ...\n'" << *sl++ << "'";
    }
    std::cout << "};" << "\n\n\n";

    auto ol = pl.output_names2str().cbegin();
    std::cout << "plant_num_outputs  = " << decltype(pl)::num_outputs << ";\n";
    std::cout << "plant_output_names = {'" << *ol++ << "'";
    while (ol != pl.output_names2str().cend()) {
        std::cout << "; ...\n'" << *ol++ << "'";
    }
    std::cout << "};" << "\n\n\n";

    auto il = pl.input_names2str().cbegin();
    std::cout << "plant_num_inputs  = " << decltype(pl)::num_inputs << ";\n";
    std::cout << "plant_input_names = {'" << *il++ << "'";
    while (il != pl.input_names2str().cend()) {
        std::cout << "; ...\n'" << *il++ << "'";
    }
    std::cout << "};" << "\n\n\n";

    auto dil = pl.aerodynamics().dataset_input_names2str().cbegin();
    std::cout << "plant_num_aerodataset_inputs  = " << std::decay_t<decltype(pl.aerodynamics())>::num_dataset_inputs << ";\n";
    std::cout << "plant_aerodataset_input_names = {'" << *dil++ << "'";
    while (dil != pl.aerodynamics().dataset_input_names2str().cend()) {
        std::cout << "; ...\n'" << *dil++ << "'";
    }
    std::cout << "};" << "\n\n\n";

    auto dcl = pl.aerodynamics().dataset_coefficient_names2str().cbegin();
    std::cout << "plant_num_aerodataset_coefficients  = " << std::decay_t<decltype(pl.aerodynamics())>::num_dataset_coefficients << ";\n";
    std::cout << "plant_aerodataset_coefficient_names = {'" << *dcl++ << "'";
    while (dcl != pl.aerodynamics().dataset_coefficient_names2str().cend()) {
        std::cout << "; ...\n'" << *dcl++ << "'";
    }
    std::cout << "};" << "\n\n\n";

    auto fml = pl.aerodynamics().force_and_moment_input_names2str().cbegin();
    std::cout << "plant_num_aeroforce_and_moment_inputs  = " << std::decay_t<decltype(pl.aerodynamics())>::num_force_and_moment_inputs << ";\n";
    std::cout << "plant_aeroforce_and_moment_input_names = {'" << *fml++ << "'";
    while (fml != pl.aerodynamics().force_and_moment_input_names2str().cend()) {
        std::cout << "; ...\n'" << *fml++ << "'";
    }
    std::cout << "};" << "\n\n\n";

    auto fmc = pl.aerodynamics().force_and_moment_coefficient_names2str().cbegin();
    std::cout << "plant_num_aeroforce_and_moment_coefficients  = " << std::decay_t<decltype(pl.aerodynamics())>::num_force_and_moment_coefficients << ";\n";
    std::cout << "plant_aeroforce_and_moment_coefficient_names = {'" << *fmc++ << "'";
    while (fmc != pl.aerodynamics().force_and_moment_coefficient_names2str().cend()) {
        std::cout << "; ...\n'" << *fmc++ << "'";
    }
    std::cout << "};" << "\n\n\n";

}
