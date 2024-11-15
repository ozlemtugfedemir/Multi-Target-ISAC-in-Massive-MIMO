function output_parameters = generate_A_B_C(input_parameters,setup_outputs)

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


%location_parameters = generate_locations(input_parameters);

positions = setup_outputs.positions;
distances = setup_outputs.distances; 
thetas = setup_outputs.thetas;
betas = setup_outputs.betas;

%plotting_locations(input_parameters,positions)


%% Channels and precoding 
%channel_vectors = cell(1,K+T);

 % for cnt_all=1:K+T
 %   if cnt_all>K %Target channels should be LOS!
 %       channel_vectors{cnt_all} = generate_channel('LOS', betas(cnt_all), thetas(cnt_all), M, noise_var);       
 %   else    %User channels 
 %       channel_vectors{cnt_all} = generate_channel('CR', betas(cnt_all), thetas(cnt_all), M, noise_var);
 %   end        
 % end
 channel_vectors = setup_outputs.channel_vectors;
 
[covariance_matrix,precoding_vector] = generate_precoding('Nullspace', channel_vectors,setup_outputs,T,K);

[feasibility_output, powers_only_UEs] = check_power_allocation_feasibility(precoding_vector, channel_vectors, K, T, SINR_th_db, P_t_max);

%[covariance_matrix_targets,precoding_vector_targets] = generate_precoding('MRT', channel_vectors(end-T+1:end) );

% for cnt_prec_t = 1:T
%     precoding_vector{end-T+cnt_prec_t} = precoding_vector_targets{cnt_prec_t};
%     covariance_matrix{end-T+cnt_prec_t} = covariance_matrix_targets{cnt_prec_t};
% end
 
 %% Target sensing 
 Array_matrices = cell(1,T);
 Array_der_matrices = cell(1,T);
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
 
 %% Finalize A_k, B_k and C_k matrices
 A = cell(1,K+T);
 B = cell(1,K+T);
 C = cell(1,K+T);
beta_min = min(beta_targets);
 for k = 1:K+T
     R_s_k = covariance_matrix{k};
     A_k = zeros(T,T);
     B_k = zeros(T,2*T);
     C_k = zeros(2*T,2*T);
     for l = 1:T
         A_matr_l = Array_matrices{l};
         A_derv_l = Array_der_matrices{l};
         alpha_l = alpha(l);
         alpha_l_norm = alpha_l/(beta_min^(1/2));

         for p = 1:T
             A_matr_p = Array_matrices{p};
             A_derv_p = Array_der_matrices{p};
             alpha_p = alpha(p);
             alpha_p_norm = alpha_p/(beta_min^(1/2));
             %% A_k part
             a_lp = trace(A_derv_p*R_s_k*(A_derv_l'));
             %a_lp = numeric_correct_complex(a_lp,1e-10);
             a_lp = real(conj(alpha_l_norm)*alpha_p_norm*a_lp);
             A_k(l,p) = a_lp;
             
             %% B_k part
             b_lp = trace(A_matr_p*R_s_k*(A_derv_l'));
             %b_lp = numeric_correct_complex(b_lp,1e-10);
             b_lp = real((conj(alpha_l_norm)*b_lp).*[1, 1i]);
             B_k(l,2*p-1:2*p) = b_lp;
             %% C_k part
             c_lp = trace(A_matr_p*R_s_k*(A_matr_l'));
             %%%%%%%%% Modify the numeric errors 
             %c_lp = numeric_correct_complex(c_lp,1e-10);
             c_lp = real(c_lp.*([1, 1i]'*[1, 1i]));
             C_k(2*l-1:2*l,2*p-1:2*p) = c_lp;
         end
     end
     A{k} = A_k;
     B{k} = B_k;
     C{k} = C_k;
 end


output_parameters.A = A;
output_parameters.B = B;
output_parameters.C = C;
output_parameters.precoding_vector = precoding_vector;
output_parameters.positions = positions;
output_parameters.channel_vectors = channel_vectors;
output_parameters.alpha = alpha;
output_parameters.beta_min = beta_min;
output_parameters.feasibility_output = feasibility_output;
output_parameters.powers_only_UEs = powers_only_UEs;
output_parameters.Array_matrices = Array_matrices;
output_parameters.Array_der_matrices = Array_der_matrices;
output_parameters.beta_targets = beta_targets;
end