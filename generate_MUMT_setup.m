function setup_output = generate_MUMT_setup(global_variables)

M = global_variables.M;
T = global_variables.T; %Number of targets
K= global_variables.K; %Number of UEs
RCS_var_db = global_variables.RCS_var_db; %RCS variance investigated in between -40 to 0. Typically -50 to 60.
number_of_setups = global_variables.number_of_setups;
angular_sep_th = 0.2;
setup_counter = 1;
noise_var_dbm = global_variables.noise_var_dbm;
noise_var = db2pow(noise_var_dbm)/1e3;

while setup_counter <= number_of_setups

    location_parameters = generate_locations(global_variables);

    positions = location_parameters.positions;
    distances = location_parameters.distances;
    thetas = location_parameters.thetas;
    betas = location_parameters.betas;


    %%% Chech the scenario and eliminate if necessary
    degree_check =  rad2deg(thetas);
    for element_i = 1:K+T
        for element_j = 1:K+T
            control_degree(element_i,element_j) = abs(degree_check(element_i) - degree_check(element_j));
        end
    end

    control_degree = control_degree + 50.*eye(K+T,K+T); % to neglect the identical elements
    minimum_angular_sep = min(min(control_degree));
    if minimum_angular_sep > angular_sep_th
        %% Channels and precoding
        channel_vectors = cell(1,K+T);

        for cnt_all=1:K+T
            if cnt_all>K %Target channels should be LOS!
                channel_vectors{cnt_all} = generate_channel('LOS', betas(cnt_all), thetas(cnt_all), M, noise_var);
            else    %User channels
                channel_vectors{cnt_all} = generate_channel('CR', betas(cnt_all), thetas(cnt_all), M, noise_var);
            end
        end


        RCS_var = db2pow(RCS_var_db);
        RCS = sqrt(RCS_var)*(-0.5502 + 0.7521i);

        % for cnt_t = 1:T
        %     RCS = sqrt(0.5*RCS_var).*(randn(1,1)+ 1i.*randn(1,1));
        % end

        setup_output.positions = positions;
        setup_output.distances = distances;
        setup_output.thetas = thetas;
        setup_output.betas = betas;
        setup_output.channel_vectors = channel_vectors;
        setup_output.RCS = RCS;

        save(['Setups_T6/separated_T', int2str(T), '_K',int2str(K),'_',int2str(setup_counter),'.mat'],"setup_output")
        setup_counter = setup_counter+1;




    else
        disp(minimum_angular_sep)
    end




end

end