function output_parameters = generate_orthogonal_input_params(input_parameters,setup_outputs)

M = input_parameters.M;
T = input_parameters.T; %Number of targets
K= input_parameters.K; %Number of UEs
N = input_parameters.N; %Number of symbols used for sensing
carrier_freq = input_parameters.carrier_freq; %Carrier frequency
noise_var_dbm = input_parameters.noise_var_dbm; %Noise variance taken from Zinat Globecom
RCS_var_db = input_parameters.RCS_var_db; %RCS variance investigated in between -40 to 0. Typically -50 to 60.
%d_max_area = input_parameters.d_max_area; %500x500m area is assumed
%BS_position = input_parameters.BS_position;
SINR_th_db = input_parameters.SINR_th_db;
%theta_tars = input_parameters.theta_tars;
%rho_tars = input_parameters.rho_tars;
P_t_max = input_parameters.P_t_max;
%rho_UEs = input_parameters.rho_UEs;
%theta_UEs = input_parameters.theta_UEs;



noise_var = db2pow(noise_var_dbm)/1e3;


%% Generate A_k, B_k and C_k matrices

positions = setup_outputs.positions;
distances = setup_outputs.distances; 
thetas = setup_outputs.thetas;
betas = setup_outputs.betas;

%plotting_locations(input_parameters,positions)


%% Channels and precoding 
% channel_vectors = cell(1,K+T);
% 
%  for cnt_all=1:K+T
%    if cnt_all>K %Target channels should be LOS!
%        channel_vectors{cnt_all} = generate_channel('LOS', betas(cnt_all), thetas(cnt_all), M, noise_var);       
%    else    %User channels 
%        channel_vectors{cnt_all} = generate_channel('CR', betas(cnt_all), thetas(cnt_all), M, noise_var);
%    end        
%  end
 
channel_vectors = setup_outputs.channel_vectors; 
[~,precoding_vector] = generate_precoding('Nullspace', channel_vectors(1:K),setup_outputs,T,K);

%[feasibility_output, powers_only_UEs] = check_orthogonal_power_allocation_feasibility(precoding_vector, channel_vectors, K, T, SINR_th_db, P_t_max);

%[covariance_matrix_targets,precoding_vector_targets] = generate_precoding('MRT', channel_vectors(end-T+1:end) );

% for cnt_prec_t = 1:T
%     precoding_vector{end-T+cnt_prec_t} = precoding_vector_targets{cnt_prec_t};
%     covariance_matrix{end-T+cnt_prec_t} = covariance_matrix_targets{cnt_prec_t};
% end
 
 %% Target sensing 
 Array_matrices = cell(1,T);
 Array_der_matrices = cell(1,T);

 RCS_var = db2pow(RCS_var_db);
 alpha = zeros(T,1);
 %alpha_tilde = zeros(2,T);
 
 for cnt_t = 1:T
    target_PL = generate_roundtrip_loss(carrier_freq, distances(cnt_t+K));
    beta_targets(cnt_t) = target_PL;
    RCS = setup_outputs.RCS;
    alpha(cnt_t) = RCS*sqrt(target_PL); %Complex Tx1
    %alpha_tilde(:,cnt_t) = [real(RCS*sqrt(target_PL));imag(RCS*sqrt(target_PL))]; %Vectorized real 2xT
    [Array_matrices{cnt_t}, Array_der_matrices{cnt_t}] = generate_A_matrices(M, [thetas(cnt_t+K),0]);
 end
 beta_min = min(beta_targets);

 
output_parameters.precoding_vector = precoding_vector;
output_parameters.positions = positions;
output_parameters.channel_vectors = channel_vectors;
output_parameters.alpha = alpha;
output_parameters.beta_min = beta_min;
%output_parameters.feasibility_output = feasibility_output;
%output_parameters.powers_only_UEs = powers_only_UEs;
output_parameters.Array_matrices = Array_matrices;
output_parameters.Array_der_matrices = Array_der_matrices;
output_parameters.beta_targets = beta_targets;
end