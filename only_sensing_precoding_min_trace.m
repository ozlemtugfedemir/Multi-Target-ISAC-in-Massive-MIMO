function optimal_output = only_sensing_precoding_min_trace(global_variables, sensing_parameters)

%% CVX

T = global_variables.T; %Number of targets
K= global_variables.K; %Number of UEs
M = global_variables.M; %Number of antennas
N = global_variables.N; %Number of antennas
noise_var_dbm = global_variables.noise_var_dbm;
noise_var = db2pow(noise_var_dbm)/1e3;

SINR_th_db = global_variables.SINR_th_db;
SINR_th = db2pow(SINR_th_db);
P_t_max = global_variables.P_t_max;
RCS_var_db = global_variables.RCS_var_db;

normalization_exp = -round((60+RCS_var_db)/10);
normalization_coeff = 10^normalization_exp;

Array_matrices = sensing_parameters.Array_matrices;
Array_der_matrices = sensing_parameters.Array_der_matrices;
beta_targets = sensing_parameters.beta_targets;
alpha = sensing_parameters.alpha;
channel_vectors = sensing_parameters.channel_vectors;

A = sensing_parameters.A;
B = sensing_parameters.B;
C = sensing_parameters.C;
precoding_vectors = sensing_parameters.precoding_vector;
channel_vectors = sensing_parameters.channel_vectors;

precoding_gain = zeros(K,K);

for cnt_ues_1 = 1:K
    precoding = precoding_vectors{cnt_ues_1};
    for cnt_ues_2 = 1:K
        channel = channel_vectors{cnt_ues_2};
        precoding_gain(cnt_ues_1,cnt_ues_2) = abs((channel)'*precoding)^2;
    end
end


cvx_begin 

variable p(K,1) nonnegative;
variable R_s(M,M) hermitian semidefinite;
variable D(T,T) symmetric;


C_sum = zeros(2*T,2*T);
B_sum = zeros(T,2*T);
A_sum = zeros(T,T);

beta_min = min(beta_targets);

C_sum = zeros(2*T,2*T);
B_sum = zeros(T,2*T);
A_sum = zeros(T,T);   

for i = 1: K
    C_sum = C_sum + p(i).*C{i};
    B_sum = B_sum + p(i).*B{i};
    A_sum = A_sum + p(i).*A{i};
end
 



    for l = 1:T
        A_matr_l = Array_matrices{l};
        A_derv_l = Array_der_matrices{l};
        alpha_l = alpha(l);
        alpha_l_norm = alpha_l/sqrt(beta_min);
        
        for z = 1:T
            A_matr_p = Array_matrices{z};
            A_derv_p = Array_der_matrices{z};
            alpha_p = alpha(z);
            alpha_p_norm = alpha_p/sqrt(beta_min);
            %% A_k part
            a_lp = trace(A_derv_p*R_s*(A_derv_l'));
            a_lp = real(conj(alpha_l_norm)*alpha_p_norm*a_lp);
            A_s(l,z) = a_lp;
            
            %% B_k part
            b_lp = trace(A_matr_p*R_s*(A_derv_l'));
            b_lp = real((conj(alpha_l_norm)*b_lp).*[1, 1i]);
            B_s(l,2*z-1:2*z) = b_lp;
            %% C_k part
            c_lp = trace(A_matr_p*R_s*(A_matr_l'));
            c_lp = real(c_lp.*([1, 1i]'*[1, 1i]));
            C_s(2*l-1:2*l,2*z-1:2*z) = c_lp;
        end
    end
    
    
C_sum = C_sum + C_s;
B_sum = B_sum + B_s;
A_sum = A_sum + A_s;
    
A_sum = 1e-4.*A_sum;
B_sum = 1e-4.*1e1.*B_sum;
C_sum = 1e-4.*1e2.*C_sum;


minimize( 1e1*trace_inv(1e+0.*D) );

subject to

%% SINR Constraints
for k = 1: K
    channel_ue = channel_vectors{k};   
    %real(channel_ue'*R_ue(:,:,k)*channel_ue)/SINR_th - real(s_int(k)) -  real(channel_ue'*R_s*channel_ue) - 1 >= 1;
    ((SINR_th+1)/SINR_th)*precoding_gain(k,k)*p(k) - precoding_gain(:,k)'*p -  real(channel_ue'*R_s*channel_ue) - 1 >= 0;
   
end

%% Sensing Constraints
[A_sum - D, B_sum; transpose(B_sum), C_sum] + 1e-3.*eye(3*T)  == semidefinite(3*T);


%% Total sum power constraint

sum(p) + trace(R_s) <= P_t_max;

cvx_end


if strcmp(cvx_status, 'Solved')
    optimal_output.p = p;
    optimal_output.R_s = R_s;
    optimal_output.D = D;
    optimal_output.status = 1;
    optimal_output.cvx_optval = cvx_optval;
    optimal_output.elapsed_time = cvx_cputime;
    
else
    optimal_output.status = 0;
    optimal_output.cvx_status = cvx_status;
    optimal_output.alpha = alpha;
    optimal_output.beta_targets = beta_targets;
end

end


