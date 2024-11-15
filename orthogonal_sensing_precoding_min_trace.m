function optimal_output = orthogonal_sensing_precoding_min_trace(global_variables, sensing_parameters)

%% CVX

T = global_variables.T; %Number of targets
K= global_variables.K; %Number of UEs
M = global_variables.M; %Number of antennas
SINR_th_db = global_variables.SINR_th_db;
SINR_th = db2pow(SINR_th_db);
P_t_max = global_variables.P_t_max;
eta = global_variables.eta;
Array_matrices = sensing_parameters.Array_matrices;
Array_der_matrices = sensing_parameters.Array_der_matrices;
beta_targets = sensing_parameters.beta_targets;
alpha = sensing_parameters.alpha;
precoding_vectors = sensing_parameters.precoding_vector;
channel_vectors = sensing_parameters.channel_vectors;


for cnt_ues_1 = 1:K
    precoding = precoding_vectors{cnt_ues_1};
    for cnt_ues_2 = 1:K
        channel = channel_vectors{cnt_ues_2};
        precoding_gain(cnt_ues_1,cnt_ues_2) = abs((channel)'*precoding)^2;
    end
end

%% Orthonal Sharing Modifications

cvx_begin 

variable R_s(M,M) hermitian semidefinite;
variable D(T,T) symmetric;
variable p(K,1) nonnegative;



beta_min = min(beta_targets);

C_sum = zeros(2*T,2*T);
B_sum = zeros(T,2*T);
A_sum = zeros(T,T);   

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
        
A_s = 1e-6.*A_s;
B_s = 1e-6.*1e1.*B_s;
C_s = 1e-6.*1e2.*C_s;


minimize( 1e2*trace_inv(1e1.*D) );

subject to


%% Sensing Constraints
[A_s - D, B_s; transpose(B_s), C_s] + 1e-3.*eye(3*T) == semidefinite(3*T);


%% Communication Constraints
for k = 1: K
((SINR_th+1)/SINR_th)*precoding_gain(k,k)*p(k) - precoding_gain(:,k)'*p - 1 >= 0;
end
        

%% Total sum power constraint

(1-eta)*trace(R_s) + eta*sum(p)  <= P_t_max;


cvx_end


if strcmp(cvx_status, 'Solved')
    optimal_output.R_s = R_s;
    optimal_output.D = D;
    optimal_output.p = p;
    optimal_output.status = 1;
    optimal_output.cvx_optval = cvx_optval;
    optimal_output.elapsed_time = cvx_cputime;
    optimal_output.cvx_status = cvx_status;
    optimal_output.alpha = alpha;
    optimal_output.beta_targets = beta_targets;
else
    optimal_output.status = 0;
    optimal_output.cvx_status = cvx_status;
    optimal_output.alpha = alpha;
    optimal_output.beta_targets = beta_targets;
end

end


