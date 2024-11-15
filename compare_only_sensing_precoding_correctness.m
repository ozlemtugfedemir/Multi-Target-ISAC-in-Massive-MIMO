function verification_output  = compare_only_sensing_precoding_correctness(optimal_output, output_parameters, global_variables)
%% Variable unbundle 
T = global_variables.T; %Number of targets
K= global_variables.K; %Number of UEs
SINR_th_db = global_variables.SINR_th_db;
SINR_th = db2pow(SINR_th_db);
p = optimal_output.p;
R_s =  optimal_output.R_s;

channel_vectors = output_parameters.channel_vectors;
precoding_vectors = output_parameters.precoding_vector;

p(p<=1e-10) = 0;


precoding_gain = zeros(K,K);

for cnt_ues_1 = 1:K
    precoding = precoding_vectors{cnt_ues_1};
    for cnt_ues_2 = 1:K
        channel = channel_vectors{cnt_ues_2};
        precoding_gain(cnt_ues_1,cnt_ues_2) = abs((channel)'*precoding)^2;
    end
end






%% Verification of SINR constraints 
for k = 1: K
    channel_ue = channel_vectors{k}; 
    SINR(k) = (precoding_gain(k,k)*p(k))/(precoding_gain(:,k)'*p - precoding_gain(k,k)*p(k) +  real(channel_ue'*R_s*channel_ue + 1));   
end



%% Sensing constraint 


verification_output.SINR = SINR;
end