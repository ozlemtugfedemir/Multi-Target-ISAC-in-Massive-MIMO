function cramer_output = calculate_crb_matr(optimal_outputs, output_parameters, global_variables)
theta_tars = global_variables.theta_tars;
K = global_variables.K;
N = global_variables.N;
T = global_variables.T;
M = global_variables.M;
alpha = output_parameters.alpha;
noise_var_dbm = global_variables.noise_var_dbm;
noise_var = db2pow(noise_var_dbm)/1e3;
A = output_parameters.A;
B = output_parameters.B;
C = output_parameters.C;
beta_min = output_parameters.beta_min;
p = optimal_outputs.p;


    A_sum = zeros(T,T);
    C_sum = zeros(2*T,2*T);
    B_sum = zeros(T,2*T);
   % B_transpose_sum = zeros(2*T,T);
        
    for i = 1: K+T
        A_sum = A_sum + p(i).*A{i};
        C_sum = C_sum + p(i).*C{i};
        B_sum = B_sum + p(i).*B{i};
       % B_transpose_sum = B_transpose_sum + p(i).*transpose(B{i});
    end
constant = 2*N/noise_var;

A_sum_unnormalized = constant*beta_min*A_sum;
B_sum_unnormalized = constant*sqrt(beta_min)*B_sum;
C_sum_unnormalized = constant*C_sum;

DOA_mat = A_sum_unnormalized - B_sum_unnormalized*inv(C_sum_unnormalized)*transpose(B_sum_unnormalized);
normalized_DOA_mat = (A_sum - B_sum*inv(C_sum)*transpose(B_sum));



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