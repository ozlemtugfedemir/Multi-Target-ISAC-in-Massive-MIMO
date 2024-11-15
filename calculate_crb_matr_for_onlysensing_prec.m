function cramer_output = calculate_crb_matr_for_onlysensing_prec(optimal_outputs, output_parameters, global_variables)
K = global_variables.K;
N = global_variables.N;
T = global_variables.T;
M = global_variables.M;
noise_var_dbm = global_variables.noise_var_dbm;
noise_var = db2pow(noise_var_dbm)/1e3;
p = optimal_outputs.p;
R_s =  optimal_outputs.R_s;

RCS_var_db = global_variables.RCS_var_db; %RCS variance investigated in between -40 to 0. Typically -50 to 60.


Array_matrices = output_parameters.Array_matrices;
Array_der_matrices = output_parameters.Array_der_matrices;
beta_targets = output_parameters.beta_targets;
alpha = output_parameters.alpha;
channel_vectors = output_parameters.channel_vectors;
A = output_parameters.A;
B = output_parameters.B;
C = output_parameters.C;


C_sum = zeros(2*T,2*T);
B_sum = zeros(T,2*T);
A_sum = zeros(T,T);

beta_min = min(beta_targets);

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
        
        for p = 1:T
            A_matr_p = Array_matrices{p};
            A_derv_p = Array_der_matrices{p};
            alpha_p = alpha(p);
            alpha_p_norm = alpha_p/sqrt(beta_min);
            %% A_k part
            a_lp = trace(A_derv_p*R_s*(A_derv_l'));
            a_lp = real(conj(alpha_l_norm)*alpha_p_norm*a_lp);
            A_s(l,p) = a_lp;
            
            %% B_k part
            b_lp = trace(A_matr_p*R_s*(A_derv_l'));
            b_lp = real((conj(alpha_l_norm)*b_lp).*[1, 1i]);
            B_s(l,2*p-1:2*p) = b_lp;
            %% C_k part
            c_lp = trace(A_matr_p*R_s*(A_matr_l'));
            c_lp = real(c_lp.*([1, 1i]'*[1, 1i]));
            C_s(2*l-1:2*l,2*p-1:2*p) = c_lp;
        end
    end
    
    
C_sum = C_sum + C_s;
B_sum = B_sum + B_s;
A_sum = A_sum + A_s;


    
constant = 2*N/noise_var;

A_sum_unnormalized = constant*beta_min*A_sum;
B_sum_unnormalized = constant*sqrt(beta_min)*B_sum;
C_sum_unnormalized = constant*C_sum;

DOA_mat = A_sum_unnormalized - B_sum_unnormalized*(C_sum_unnormalized\transpose(B_sum_unnormalized));
normalized_DOA_mat = (A_sum - B_sum*(C_sum\transpose(B_sum)));



normalized_CRB_mat = inv(normalized_DOA_mat);
CRB_mat = inv(DOA_mat);


for cnt_tars = 1:T
SNR = N*abs(alpha(cnt_tars))^2/noise_var;
SNR_tars(cnt_tars) = SNR;
end


cramer_output.CRB = diag(CRB_mat);
cramer_output.normalized_CRB = diag(normalized_CRB_mat);
cramer_output.SNR_tars = SNR_tars;



end