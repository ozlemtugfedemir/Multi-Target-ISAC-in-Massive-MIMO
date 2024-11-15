function optimal_output = precoding_min_trace(global_variables, sensing_parameters)

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



Array_matrices = sensing_parameters.Array_matrices;
Array_der_matrices = sensing_parameters.Array_der_matrices;
beta_targets = sensing_parameters.beta_targets;
alpha = sensing_parameters.alpha;
channel_vectors = sensing_parameters.channel_vectors;


cvx_begin 

variable R_ue(M,M,K) hermitian semidefinite;
variable R_s(M,M) hermitian semidefinite;
variable D(T,T) symmetric;


C_sum = zeros(2*T,2*T);
B_sum = zeros(T,2*T);
A_sum = zeros(T,T);

R_s_k = zeros(M,M);
beta_min = min(beta_targets);

normalization_exp = -round((60+RCS_var_db)/10);
normalization_coeff = 10^normalization_exp;

for cnt_k = 1:K+1
    if cnt_k <= K
        R_s_k = R_ue(:,:,cnt_k);
    else
        R_s_k = R_s;
    end
%     A_k = zeros(T,T);
%     B_k = zeros(T,2*T);
%     C_k = zeros(2*T,2*T);
    for l = 1:T
        A_matr_l = Array_matrices{l};
        A_derv_l = Array_der_matrices{l};
        alpha_l = alpha(l);
        alpha_l_norm = alpha_l/sqrt(beta_min);
        
        for p = 1:T
            A_matr_p = Array_matrices{p};
            A_derv_p = Array_der_matrices{p};
            alpha_p = alpha(p);
            alpha_p_norm = alpha_p/sqrt(beta_min);
            %% A_k part
            a_lp = trace(A_derv_p*R_s_k*(A_derv_l'));
            a_lp = real(conj(alpha_l_norm)*alpha_p_norm*a_lp);
            A_k(l,p) = a_lp;
            
            %% B_k part
            b_lp = trace(A_matr_p*R_s_k*(A_derv_l'));
            b_lp = real((conj(alpha_l_norm)*b_lp).*[1, 1i]);
            B_k(l,2*p-1:2*p) = b_lp;
            %% C_k part
            c_lp = trace(A_matr_p*R_s_k*(A_matr_l'));
            c_lp = real(c_lp.*([1, 1i]'*[1, 1i]));
            C_k(2*l-1:2*l,2*p-1:2*p) = c_lp;
        end
    end
    A_sum = A_k + A_sum;
    B_sum = B_k + B_sum;
    C_sum = C_k + C_sum;
end

   A_sum = 1e-7.*A_sum;
   B_sum = 1e1.*1e-7.*B_sum;
   C_sum = 1e2*1e-7.*C_sum;

%A_sum = constant^-1.*A_sum;
%B_sum = constant^-1.*B_sum;
%C_sum = constant^-1.*C_sum;

%M_matr = 1e-2.*[A_sum, B_sum; transpose(B_sum), C_sum];

minimize( 1e1.*trace_inv(1e3.*D) );

subject to

%% SINR Constraints
for k = 1: K
    
    channel_ue = channel_vectors{k};   
    for cnt_interf = 1:K
        sum_intf_gain = 0;
        if cnt_interf == k
            intf_gain = 0;
        else
            intf_gain = channel_ue'*R_ue(:,:,cnt_interf)*channel_ue;
        end
        sum_intf_gain = intf_gain+sum_intf_gain;
    end
    s_int(k) = sum_intf_gain;
    real(channel_ue'*R_ue(:,:,k)*channel_ue)/SINR_th - real(s_int(k)) -  real(channel_ue'*R_s*channel_ue) - 1 >= 1;
   
end

%% Sensing Constraints
[A_sum - D, B_sum; transpose(B_sum), C_sum] + 1e-3.*eye(3*T) == semidefinite(3*T);



%% Total sum power constraint

%tr_ue = zeros(K,1);
for cnt_u = 1:K
   tr_ue(cnt_u) = trace(R_ue(:,:,cnt_u));
end

sum(tr_ue) + trace(R_s) <= P_t_max;
 






cvx_end


if strcmp(cvx_status, 'Solved')
    optimal_output.R_ue = R_ue;
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


