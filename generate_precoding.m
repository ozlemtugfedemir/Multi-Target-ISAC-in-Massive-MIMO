function [covariance_matrix,precoding_vector] = generate_precoding(precoding_type,channel_vectors,setup_outputs,T,K)
%%%%% Provides the unit-norm precoding vector. Options for now: MRT and RZF. 
%%%%% channel_vectors is a cell-array holding the channel coefficient of
%%%%% UEs
number_of_UEs = length(channel_vectors);
number_of_antennas = length(channel_vectors{1});

precoding_vector = cell(1,number_of_UEs);
covariance_matrix = cell(1,number_of_UEs);


    if strcmp(precoding_type,'MRT')
        
        for cnt = 1:number_of_UEs 
            UE_channel = channel_vectors{cnt};
            precoding_vector_for_UE = (UE_channel)./(norm(UE_channel));
            precoding_vector{cnt} = precoding_vector_for_UE;
            covariance_matrix{cnt} = precoding_vector_for_UE*precoding_vector_for_UE';
        end
            
    elseif strcmp(precoding_type,'RZF')
            H_matrix = zeros(number_of_antennas, number_of_UEs);
            lambda = 0.01;
            %% Unbundle the channel matrix
            for cnt = 1:number_of_UEs
                H_matrix(:,cnt) = channel_vectors{cnt};
            end    
            precoding_matrix = H_matrix*inv(H_matrix'*H_matrix+lambda.*eye(number_of_UEs));
            
            %% Bundle the precoding 
            for cnt_2 = 1:number_of_UEs
               precoding_vector_for_UE  = precoding_matrix(:,cnt_2);
               precoding_vector_for_UE = precoding_vector_for_UE./(norm(precoding_vector_for_UE));
               precoding_vector{cnt_2}  = precoding_vector_for_UE;
               covariance_matrix{cnt_2} = precoding_vector_for_UE*precoding_vector_for_UE';
            end  
    elseif    strcmp(precoding_type,'Nullspace')
            H_matrix = zeros(number_of_antennas, K);
            %% Unbundle the channel matrix
            for cnt = 1:K
                H_matrix(:,cnt) = channel_vectors{cnt};
            end
            %% R-ZF for UEs
            lambda = 0.01;
            precoding_matrix = H_matrix*inv(H_matrix'*H_matrix+lambda.*eye(K));
           %% Bundle the precoding 
            for cnt_2 = 1:K
               precoding_vector_for_UE  = precoding_matrix(:,cnt_2);
               precoding_vector_for_UE = precoding_vector_for_UE./(norm(precoding_vector_for_UE));
               precoding_vector{cnt_2}  = precoding_vector_for_UE;
               covariance_matrix{cnt_2} = precoding_vector_for_UE*precoding_vector_for_UE';
            end              
            
            
            %% Nullspace projection for targets
            [U,~,~] = svd(H_matrix);
            U_prod = U*U';
            prod_matrix = eye(size(U_prod))-U_prod;
            thetas = setup_outputs.thetas;
            for cnt_t =1:T
            theta_tar = thetas(K+cnt_t);
            [~,arr_res_t] = generate_array_response('ULA', number_of_antennas, [theta_tar,0]);
            precoding_vector_for_target = (prod_matrix*arr_res_t)./(norm((prod_matrix*arr_res_t)));
            precoding_vector{K+cnt_t}  = precoding_vector_for_target;
            covariance_matrix{K+cnt_t} = precoding_vector_for_target*precoding_vector_for_target';
            end
    end
end
