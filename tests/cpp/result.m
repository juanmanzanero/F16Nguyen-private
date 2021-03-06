plant_num_states  = 21;
plant_state_names = {'u_mps'; ...
'v_mps'; ...
'w_mps'; ...
'p_radps'; ...
'q_radps'; ...
'r_radps'; ...
'qw_body2Earth'; ...
'qx_body2Earth'; ...
'qy_body2Earth'; ...
'qz_body2Earth'; ...
'xEarth_m'; ...
'yEarth_m'; ...
'zEarth_m'; ...
'dh_deg'; ...
'dlef_deg'; ...
'dsb_deg'; ...
'da_deg'; ...
'dr_deg'; ...
'thrustvectoring_longitude_angle_deg'; ...
'thrustvectoring_latitude_angle_deg'; ...
'P3_percent'};


plant_num_outputs  = 104;
plant_output_names = {'dynamic_pressure_Pa'; ...
'Mach'; ...
'VTASx_mps'; ...
'VTASy_mps'; ...
'VTASz_mps'; ...
'TAS_mps'; ...
'aoa_deg'; ...
'aos_deg'; ...
'true_heading_angle_deg'; ...
'flight_path_angle_deg'; ...
'flight_roll_angle_deg'; ...
'pw_radps'; ...
'qw_radps'; ...
'rw_radps'; ...
'paero_radps'; ...
'qaero_radps'; ...
'raero_radps'; ...
'air_density_kgpm3'; ...
'sounspeed_mps'; ...
'air_kinematic_viscosity_m2ps'; ...
'air_temperature_K'; ...
'air_pressure_Pa'; ...
'Faerox_N'; ...
'Faeroy_N'; ...
'Faeroz_N'; ...
'Maerox_Nm'; ...
'Maeroy_Nm'; ...
'Maeroz_Nm'; ...
'CXtot'; ...
'CYtot'; ...
'CZtot'; ...
'CLLtot'; ...
'CMtot'; ...
'CNtot'; ...
'engine_thrust_N'; ...
'engine_angular_momentum_kgm2ps'; ...
'P1_percent'; ...
'P2_percent'; ...
'invtau_1ps'; ...
'Tidle_N'; ...
'Tmil_N'; ...
'Tmax_N'; ...
'VGSx_Earthaxes_mps'; ...
'VGSy_Earthaxes_mps'; ...
'VGSz_Earthaxes_mps'; ...
'DOWNx_bodyaxes'; ...
'DOWNy_bodyaxes'; ...
'DOWNz_bodyaxes'; ...
'yawdot_degps'; ...
'pitchdot_degps'; ...
'rolldot_degps'; ...
'yaw_deg'; ...
'pitch_deg'; ...
'roll_deg'; ...
'GS_mps'; ...
'ground_aoa_deg'; ...
'ground_aos_deg'; ...
'ground_true_heading_angle_deg'; ...
'ground_flight_path_angle_deg'; ...
'ground_flight_roll_angle_deg'; ...
'Ftotx_N'; ...
'Ftoty_N'; ...
'Ftotz_N'; ...
'Mtotx_N'; ...
'Mtoty_N'; ...
'Mtotz_N'; ...
'Nx_g'; ...
'Ny_g'; ...
'Nz_g'; ...
'ax_mps2'; ...
'ay_mps2'; ...
'az_mps2'; ...
'udot_mps2'; ...
'vdot_mps2'; ...
'wdot_mps2'; ...
'pdot_radps2'; ...
'qdot_radps2'; ...
'rdot_radps2'; ...
'GSdot_mps'; ...
'ground_aoadot_degps'; ...
'ground_aosdot_degps'; ...
'ground_true_heading_angledot_degps'; ...
'ground_flight_path_angledot_degps'; ...
'ground_flight_roll_angledot_degps'; ...
'ground_turn_radius_m'; ...
'dVTASdtx_mps2'; ...
'dVTASdty_mps2'; ...
'dVTASdtz_mps2'; ...
'VTASxdot_mps2'; ...
'VTASydot_mps2'; ...
'VTASzdot_mps2'; ...
'TASdot_mps2'; ...
'aoadot_degps'; ...
'aosdot_degps'; ...
'Nwx_g'; ...
'Nwy_g'; ...
'Nwz_g'; ...
'p_windREarth_radps'; ...
'q_windREarth_radps'; ...
'r_windREarth_radps'; ...
'true_heading_angledot_degps'; ...
'flight_path_angledot_degps'; ...
'flight_roll_angledot_degps'; ...
'air_turn_radius_m'};


plant_num_inputs  = 9;
plant_input_names = {'dh_dmd_deg'; ...
'dlef_dmd_deg'; ...
'dsb_dmd_deg'; ...
'da_dmd_deg'; ...
'dr_dmd_deg'; ...
'thrustvectoring_longitude_angle_dmd_deg'; ...
'thrustvectoring_latitude_angle_dmd_deg'; ...
'throttle_percent'; ...
'aoadot_degps'};


plant_num_aerodataset_inputs  = 3;
plant_aerodataset_input_names = {'aoa_deg'; ...
'aos_deg'; ...
'dh_deg'};


plant_num_aerodataset_coefficients  = 55;
plant_aerodataset_coefficient_names = {'CX'; ...
'CXdh0'; ...
'CXlef'; ...
'CXq'; ...
'DCXsb'; ...
'DCXqlef'; ...
'CY'; ...
'CYda20'; ...
'CYdr30'; ...
'CYlef'; ...
'CYda20lef'; ...
'CYp'; ...
'CYr'; ...
'DCYplef'; ...
'DCYrlef'; ...
'CZ'; ...
'CZdh0'; ...
'CZlef'; ...
'CZq'; ...
'DCZsb'; ...
'DCZqlef'; ...
'CLL'; ...
'CLLdh0'; ...
'CLLda20'; ...
'CLLdr30'; ...
'CLLlef'; ...
'CLLda20lef'; ...
'CLLp'; ...
'CLLr'; ...
'DCLLbeta'; ...
'DCLLplef'; ...
'DCLLrlef'; ...
'CM'; ...
'CMdh0'; ...
'CMlef'; ...
'DCMds'; ...
'CMq'; ...
'DCM'; ...
'DCMsb'; ...
'DCMqlef'; ...
'CN'; ...
'CNdh0'; ...
'CNda20'; ...
'CNdr30'; ...
'CNlef'; ...
'CNda20lef'; ...
'CNp'; ...
'CNr'; ...
'DCNbeta'; ...
'DCNda'; ...
'DCNplef'; ...
'DCNrlef'; ...
'etadh'; ...
'CLalphadot'; ...
'CMalphadot'};


plant_num_aeroforce_and_moment_inputs  = 15;
plant_aeroforce_and_moment_input_names = {'TAS_mps'; ...
'aoa_deg'; ...
'aos_deg'; ...
'dh_deg'; ...
'dlef_deg'; ...
'dsb_deg'; ...
'paero_radps'; ...
'qaero_radps'; ...
'raero_radps'; ...
'da_deg'; ...
'dr_deg'; ...
'aoadot_radps'; ...
'xcg_perMAC'; ...
'MAC_m'; ...
'wingspan_m'};


plant_num_aeroforce_and_moment_coefficients  = 6;
plant_aeroforce_and_moment_coefficient_names = {'CXtot'; ...
'CYtot'; ...
'CZtot'; ...
'CLLtot'; ...
'CMtot'; ...
'CNtot'};


