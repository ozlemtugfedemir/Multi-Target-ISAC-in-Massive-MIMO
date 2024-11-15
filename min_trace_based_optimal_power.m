function optimal_output = min_trace_based_optimal_power(global_variables, output_parameters)

%% CVX 

T = global_variables.T; %Number of targets
K= global_variables.K; %Number of UEs
SINR_th_db = global_variables.SINR_th_db;
SINR_th = db2pow(SINR_th_db);
P_t_max = global_variables.P_t_max;
RCS_var_db = global_variables.RCS_var_db;

A = output_parameters.A;
B = output_parameters.B;
C = output_parameters.C;
precoding_vectors = output_parameters.precoding_vector;
channel_vectors = output_parameters.channel_vectors;

normalization_exp = -round((60+RCS_var_db)/10);
normalization_coeff = 10^normalization_exp;

precoding_gain = zeros(K+T,K);

for cnt_ues_1 = 1:K+T
    precoding = precoding_vectors{cnt_ues_1};
    for cnt_ues_2 = 1:K
        channel = channel_vectors{cnt_ues_2};
        precoding_gain(cnt_ues_1,cnt_ues_2) = abs((channel)'*precoding)^2;
    end
end


cvx_begin 

    variable p(K+T,1) nonnegative;
    variable D(T,T) symmetric; 
    C_sum = zeros(2*T,2*T);
    B_sum = zeros(T,2*T);
    A_sum = zeros(T,T);    
    for i = 1: K+T
        C_sum = C_sum + p(i).*C{i};
        B_sum = B_sum + p(i).*B{i};
        A_sum = A_sum + p(i).*A{i};
    end
 
   A_sum = normalization_coeff.*A_sum;
   B_sum = 1e1.*normalization_coeff.*B_sum;
   C_sum = 1e2*normalization_coeff.*C_sum;
    
  %   A_sum = 2e-1.*1e-3.*A_sum;
 %  B_sum = 1e-3.*B_sum;
 %  C_sum = 2e1.*1e-3.*C_sum;
    
    minimize( 1e1.*trace_inv( D ) );
    
    subject to
    
        %% SINR Constraints 
        for k = 1: K
        ((SINR_th+1)/SINR_th)*precoding_gain(k,k)*p(k) - precoding_gain(:,k)'*p - 1 >= 0;
        end
        
        %% Sensing Constraints 
        [(A_sum) - D, B_sum; transpose(B_sum), C_sum] +1e-3.*eye(3*T,3*T) == semidefinite(3*T);
    
        
        %% Total sum power constraint 
        
        sum(p) <= P_t_max;
       
      
       



cvx_end


if strcmp(cvx_status, 'Solved')
    optimal_output.p = p;
    optimal_output.D = D;
    optimal_output.status = 1;
    optimal_output.cvx_optval = cvx_optval;
    optimal_output.elapsed_time = cvx_cputime;
else
    optimal_output.status = 0;   
    optimal_output.cvx_status = cvx_status;
end

end

