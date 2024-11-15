clc; clear; close all;



%% Global variables
global_variables.M = 32; %Number of BS antennas

SINR_th_db_vec = [10];
T_vec = [7];
K_vec = [8];
global_variables.N = 100; %Number of symbols used for sensing
global_variables.BW= 20e6; %Bandwidth taken from Zinat Globecom
global_variables.carrier_freq = 1.9e9; %Carrier frequency
global_variables.noise_var_dbm = -94; %Noise variance taken from Zinat Globecom
global_variables.P_t_max = 10; %in watts taken from Zinat Globecom
global_variables.RCS_var_db = 0;

global_variables.d_max_area = 500; %500x500m area is assumed
global_variables.BS_position = [0,0];
%global_variables.rho_UEs = 250;
%global_variables.theta_UEs = [linspace(pi/6,pi/3,global_variables.K)];


global_variables.rho_tars = 150;
global_variables.RCS_var_db = 0; %RCS variance investigated in between -40 to 0. Typically -50 to 60.
global_variables.UE_position_scenario = 'random';
global_variables.target_position_scenario = 'random';

%global_variables.UE_center = [90,170];
global_variables.UE_center = [-100,90];
global_variables.UE_range = 75;
global_variables.Target_center = [-100,90];
global_variables.Target_range = 75;

%global_variables.SINR_th_db = 10;
global_variables.RCS_var_db = 0;



for cnt_case = 1:length(K_vec)
        rng(1234);
        T = T_vec(1);
        K = K_vec(1);
        SINR_th_db = SINR_th_db_vec(cnt_case);
        scenario_name = ['M32_precoding_T',int2str(T),'_K', int2str(K),'_SINRt',int2str(SINR_th_db)];
        sensing_optimal_precoding_simulation(global_variables,T,K,SINR_th_db,scenario_name);

end



    