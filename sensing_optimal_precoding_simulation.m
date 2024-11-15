function sensing_optimal_precoding_simulation(global_variables,T,K,SINR_th_db,scenario_name)

global_variables.T = T; %Number of targets
global_variables.K = K; %Number of UEs
global_variables.theta_UEs = [linspace(-pi/6,-pi/8,global_variables.K)];
global_variables.theta_tars = [linspace(-pi/3,pi/3,global_variables.T)];
global_variables.SINR_th_db = SINR_th_db;

%target_degrees = (-((T-1)*degree_intv)/2:degree_intv:((T-1)*degree_intv)/2);
%global_variables.theta_tars = [deg2rad(target_degrees)];



rng(1234);
global_variables.scenario_name = scenario_name;

num_sim = 1;
tot_sim = 1;
fail_sim = 1;
status_fail = [];
alpha_fail = [];
positions_fail = [];

while tot_sim<2
    
    output_parameters = generate_sensing_parameters(global_variables);
    optimal_outputs = precoding_min_trace(global_variables, output_parameters);
    
    if optimal_outputs.status == 1
        
        verification_output  = compare_precoding_correctness(optimal_outputs, output_parameters, global_variables);
        cramer_output = calculate_crb_matr_for_prec(optimal_outputs, output_parameters, global_variables);
        %%%%%%%% if verification output.e_99 neq 1 
        %%%%%%%% rank reduction algorithm add 
        
        
        
        
        SINRS = verification_output.SINR;
        e_9 = verification_output.e_cond_9;
        e_99 = verification_output.e_cond_99;
        e_999 = verification_output.e_cond_999;
        positions = output_parameters.positions;
        CRB_vec = cramer_output.CRB;
        normalized_CRB_vec = cramer_output.normalized_CRB;
        SNR_tars = cramer_output.SNR_tars;
        objectives = optimal_outputs.cvx_optval;
        R_ue = optimal_outputs.R_ue;
        R_s = optimal_outputs.R_s;
        D = optimal_outputs.D;
        elapsed_time = optimal_outputs.elapsed_time;
        num_sim = num_sim+1
        sim_status = 1;
        save([scenario_name,'_',int2str(tot_sim)],'SINRS', 'e_9', 'e_99', 'e_999', 'positions',   'CRB_vec','SNR_tars', 'objectives', 'R_ue' , 'R_s', 'global_variables','sim_status', 'normalized_CRB_vec', 'elapsed_time', 'D');

    else
        status_fail = optimal_outputs.cvx_status;
        alpha_fail = optimal_outputs.alpha;
        positions_fail = output_parameters.positions;
      %  SINRS(tot_sim,:) = NaN.*ones(1,global_variables.K);
       % positions(tot_sim,:,:) = NaN.*ones(global_variables.K+global_variables.T,2);
      %  CRB_vec(tot_sim,:) = NaN.*ones(1,global_variables.T);
       % normalized_CRB_vec(tot_sim,:) = NaN.*ones(1,global_variables.T);
       % SNR_tars(tot_sim,:) = NaN.*ones(1,global_variables.T);
       % objectives(tot_sim) = NaN;
       % R_ue(:,:,:,tot_sim) = NaN.*ones(global_variables.M,global_variables.M,global_variables.K);
       % R_s(:,:,tot_sim) = NaN.*ones(global_variables.M,global_variables.M);
       % elapsed_time(tot_sim) = NaN;
        sim_status = 0;
        fail_sim = fail_sim+1;
       % D(:,:,tot_sim) = NaN.*ones(global_variables.T,global_variables.T);
        save([scenario_name,'_',int2str(tot_sim)],'global_variables', 'status_fail', 'alpha_fail','positions_fail','sim_status');

    end
    tot_sim = tot_sim+1;
    clearvars cramer_output verification_output optimal_outputs output_parameters normalized_cramer_output
end



end