function [feasibility_output, powers_only_UEs] = check_power_allocation_feasibility(precoding_vectors, channel_vectors, K, T, SINR_th_db, P_t_max)

precoding_gain = zeros(K+T,K);

for cnt_ues_1 = 1:K+T
    precoding = precoding_vectors{cnt_ues_1};
    for cnt_ues_2 = 1:K
        channel = channel_vectors{cnt_ues_2};
        precoding_gain(cnt_ues_1,cnt_ues_2) = abs((channel)'*precoding)^2;
    end
end

SINR_th = db2pow(SINR_th_db);

cvx_begin quiet

    variable p(K+T,1) nonnegative;
    
    minimize( 1 );
    
    subject to
    
        %% SINR Constraints 
        for k = 1: K
        ((SINR_th+1)/SINR_th)*precoding_gain(k,k)*p(k) - precoding_gain(:,k)'*p - 1 >= 0;
        end
        
        %% Total sum power constraint 
        
        sum(p) <= P_t_max;
    
cvx_end


if strcmp(cvx_status, 'Solved')
    for k = 1: K
        powers_only_UEs = p;
    end
    feasibility_output = 1
else
    feasibility_output = 0
    powers_only_UEs = NaN.*ones(K+T,1);
end





end