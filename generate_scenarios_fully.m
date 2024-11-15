clc; clear; close all;

%%%%%%% Generate Scenario and Save all Settings 
global_variables.M = 32; %Number of BS antennas
SINR_th_db_vec = [10];
global_variables.T = 8;
global_variables.K = 16;
global_variables.N = 100; %Number of symbols used for sensing
global_variables.BW= 20e6; %Bandwidth taken from Zinat Globecom
global_variables.carrier_freq = 1.9e9; %Carrier frequency
global_variables.noise_var_dbm = -94; %Noise variance taken from Zinat Globecom
global_variables.P_t_max = 10; %in watts taken from Zinat Globecom

global_variables.d_max_area = 500; %500x500m area is assumed
global_variables.BS_position = [0,0];
%global_variables.rho_UEs = 250;
%global_variables.theta_UEs = [linspace(pi/6,pi/3,global_variables.K)];


global_variables.rho_tars = 150;
global_variables.UE_position_scenario = 'random';
global_variables.target_position_scenario = 'random';

global_variables.UE_center = [240,-100]; %Sep located
%global_variables.UE_center = [170,30]; %Co-located
global_variables.UE_range = 90;
global_variables.Target_center = [170,30];
global_variables.Target_range = 75;

%global_variables.SINR_th_db = 10;
global_variables.RCS_var_db = -20;
global_variables.number_of_setups = 2000;

%%%% Outputs only locations, channel vectors and RCS values! 


%%% Generates and saves all setups at once. 
generate_MUMT_setup(global_variables);
