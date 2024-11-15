function sensing_onlytarget_precoding_simulation(global_variables,T,K,SINR_th_db,scenario_name)

global_variables.T = T; %Number of targets
global_variables.K = K; %Number of UEs
%global_variables.theta_UEs = [linspace(-pi/6,-pi/8,global_variables.K)];
%global_variables.theta_tars = [linspace(-pi/3,pi/3,global_variables.T)];
global_variables.SINR_th_db = SINR_th_db;
scenario_type = global_variables.scenario_type; 
%target_degrees = (-((T-1)*degree_intv)/2:degree_intv:((T-1)*degree_intv)/2);
%global_variables.theta_tars = [deg2rad(target_degrees)];



%rng(1234);
global_variables.scenario_name = scenario_name;

num_sim = 1;
tot_sim = 1;
fail_sim = 1;
status_fail = [];
alpha_fail = [];
positions_fail = [];


while tot_sim<500
    load(['Setups_T6/',scenario_type,'_T', int2str(T), '_K',int2str(K),'_',int2str(tot_sim),'.mat'])
    output_parameters = generate_A_B_C(global_variables,setup_output);
    optimal_outputs = only_sensing_precoding_min_trace(global_variables, output_parameters);
    
    if optimal_outputs.status == 1
        verification_output  = compare_only_sensing_precoding_correctness(optimal_outputs, output_parameters, global_variables);
        cramer_output = calculate_crb_matr_for_onlysensing_prec(optimal_outputs, output_parameters, global_variables);
        
        SINRS = verification_output.SINR;
        positions = output_parameters.positions;
        CRB_vec = cramer_output.CRB;
        normalized_CRB_vec = cramer_output.normalized_CRB;
        SNR_tars = cramer_output.SNR_tars;
        objectives = optimal_outputs.cvx_optval;
        p = optimal_outputs.p;
        R_s = optimal_outputs.R_s;
        D = optimal_outputs.D;
        elapsed_time(tot_sim) = optimal_outputs.elapsed_time;
        num_sim = num_sim+1
        sim_status = 1;
        save([scenario_name,'_',int2str(tot_sim)],'SINRS', 'positions',   'CRB_vec','SNR_tars', 'objectives', 'p' , 'R_s', 'global_variables','sim_status', 'normalized_CRB_vec', 'elapsed_time', 'D');
        
    else
        status_fail{fail_sim} = optimal_outputs.cvx_status;
        alpha_fail(:,fail_sim) = optimal_outputs.alpha;
        positions_fail(fail_sim,:,:) = output_parameters.positions;
        sim_status = 0;
        fail_sim = fail_sim+1;
        save([scenario_name,'_',int2str(tot_sim)],'status_fail', 'alpha_fail',   'positions_fail','sim_status');

    end
    tot_sim = tot_sim+1;
    clearvars cramer_output verification_output optimal_outputs output_parameters normalized_cramer_output
end


end