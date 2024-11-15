function channel_coefficient = generate_channel(channel_type, beta_db, theta, number_of_antennas, noise_var)
%%%% channel_type = 'CR (Correlated Rayleigh)', 'LOS'
%%% beta = input from path-loss 
switch channel_type
    case  'CR'
        shadowing = 4*randn(1,1);
        beta_new_db = beta_db + shadowing;
        beta_new = db2pow(beta_new_db);
        ASDdeg = 10;
        R = functionRlocalscattering(number_of_antennas,theta,ASDdeg);
        new_R = beta_new.*R;
        channel_coefficient = sqrt(new_R)*((1/sqrt(2))*(randn(number_of_antennas,1)+ 1i.*randn(number_of_antennas,1)) );  
        
    case  'LOS'
        beta = db2pow(beta_db);
        [~,array_response] = generate_array_response('ULA', number_of_antennas, [theta,0]);
        channel_coefficient = sqrt(beta).*array_response;
    otherwise    
      warning('Unexpected channel type!')   
end

channel_coefficient = channel_coefficient./(sqrt(noise_var));
end