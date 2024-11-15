function  sensing_optimal_power_simulation(global_variables,T,K,SINR_th_db,scenario_name)

global_variables.T = T; %Number of targets
global_variables.K = K; %Number of UEs
global_variables.theta_UEs = [linspace(-pi/6,-pi/8,global_variables.K)];
global_variables.theta_tars = [linspace(-pi/3,pi/3,global_variables.T)];
global_variables.SINR_th_db = SINR_th_db;

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

        while tot_sim<200
            
            output_parameters = generate_A_B_C(global_variables); %2N/noise variance effect is cancelled out!!    
            optimal_outputs =  min_trace_based_optimal_power(global_variables, output_parameters);
            
            if optimal_outputs.status == 1
                
                verification_output  = compare_correctness(optimal_outputs, output_parameters, global_variables);
                cramer_output = calculate_crb_matr(optimal_outputs, output_parameters, global_variables);
                SINRS = verification_output.SINR;
                positions = output_parameters.positions;
                CRB_vec = cramer_output.CRB;
                normalized_CRB_vec = cramer_output.normalized_CRB;
                SNR_tars= cramer_output.SNR_tars;
                objectives = optimal_outputs.cvx_optval;
                p = optimal_outputs.p;
                D = optimal_outputs.D;
                elapsed_time(tot_sim) = optimal_outputs.elapsed_time;
                num_sim = num_sim+1
                sim_status = 1;
                save([scenario_name,'_',int2str(tot_sim)],'SINRS', 'positions',   'CRB_vec','SNR_tars', 'objectives', 'p', 'global_variables','sim_status', 'normalized_CRB_vec', 'elapsed_time', 'D');
                
            else

                status_fail = optimal_outputs.cvx_status;
                alpha_fail = output_parameters.alpha;
                positions_fail = output_parameters.positions;
                sim_status = 0;
                
                if output_parameters.feasibility_output == 1
                    optimal_outputs.p = output_parameters.powers_only_UEs;
                    cramer_output = calculate_crb_matr(optimal_outputs, output_parameters, global_variables);
                    CRB_vec = cramer_output.CRB;
                    normalized_CRB_vec = cramer_output.normalized_CRB;
                    p = optimal_outputs.p;    
                    save([scenario_name,'_',int2str(tot_sim)],'CRB_vec', 'normalized_CRB_vec', 'p', 'positions_fail','status_fail', 'alpha_fail', 'sim_status');
                else
                    save([scenario_name,'_',int2str(tot_sim)], 'positions_fail','status_fail', 'alpha_fail', 'sim_status');
                end                
                fail_sim = fail_sim+1;

            end
            tot_sim = tot_sim+1;
            clearvars cramer_output verification_output optimal_outputs output_parameters 
        end 
        

end