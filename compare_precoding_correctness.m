function verification_output  = compare_precoding_correctness(optimal_output, output_parameters, global_variables)
%% Variable unbundle 
T = global_variables.T; %Number of targets
K= global_variables.K; %Number of UEs
SINR_th_db = global_variables.SINR_th_db;
SINR_th = db2pow(SINR_th_db);
R_ue = optimal_output.R_ue;
R_s =  optimal_output.R_s;

channel_vectors = output_parameters.channel_vectors;





e_thr_1 = 1-1e-3;
e_thr_2 = 1-1e-2;
e_thr_3 = 1-1e-1;

for cnt_ues_1 = 1:K
    e = eig(R_ue(:,:,cnt_ues_1));
    e_rat = max(e)/sum(e);
    e_cond_99(cnt_ues_1) = e_rat>=e_thr_2;
    e_cond_999(cnt_ues_1) = e_rat>=e_thr_1;
    e_cond_9(cnt_ues_1) = e_rat>=e_thr_3;
end

e_t = eig(R_s);
e_rat_t = max(e_t)/sum(e_t);
e_cond_99(K+1) = e_rat_t>=e_thr_2;
e_cond_999(K+1) = e_rat_t>=e_thr_1;
e_cond_9(K+1) = e_rat_t>=e_thr_3;



%% Verification of SINR constraints 
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
    SINR(k) = real(channel_ue'*R_ue(:,:,k)*channel_ue) / (real(sum_intf_gain) +  real(channel_ue'*R_s*channel_ue + 1));
   
end



%% Sensing constraint 


% for t = 1: T    
%  sensing_const(t) = s(t) - transpose(identity(:,t))*B_sum*inv(C_sum)*transpose(B_sum)*identity(:,t) >= -1e-5;     
% end
verification_output.e_cond_9 = e_cond_9;
verification_output.e_cond_99 = e_cond_99;
verification_output.e_cond_999 = e_cond_999;
verification_output.SINR = SINR;
end