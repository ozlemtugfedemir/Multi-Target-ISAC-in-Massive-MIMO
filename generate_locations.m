function location_parameters = generate_locations(input_parameters)
M = input_parameters.M;
T = input_parameters.T; %Number of targets
K= input_parameters.K; %Number of UEs
d_max_area = input_parameters.d_max_area; %500x500m area is assumed
BS_position = input_parameters.BS_position;





UE_position_scenario = input_parameters.UE_position_scenario;
target_position_scenario = input_parameters.target_position_scenario;

positions = zeros(K+T,2); 
distances = zeros(K+T,1); 
thetas = zeros(K+T,1);
betas = zeros(K+T,1);



if strcmp(UE_position_scenario,'random') && strcmp(target_position_scenario,'random')
    UE_center = input_parameters.UE_center;
    UE_range = input_parameters.UE_range;
    Target_center = input_parameters.Target_center;
    Target_range = input_parameters.Target_range;
    
        %% Generate UE positions    
    for cnt_ue = 1:K
        UE_position = (UE_range/2)*(2*rand(1,2)-1)+UE_center; %x-y
        UE_position = UE_position + BS_position;
        positions(cnt_ue,:) = UE_position;
        distance_UE = norm(UE_position-BS_position);
        distances(cnt_ue) = distance_UE;
        theta_UE = angle((UE_position-BS_position)*[1;1i]);
        thetas(cnt_ue) = theta_UE;
        beta_db_UE = generate_path_loss(distance_UE);
        betas(cnt_ue) = beta_db_UE;
    end

    for cnt_tar = 1:T
        tar_position = (Target_range/2).*(2*rand(1,2)-1)+Target_center; %x-y
        tar_position = tar_position + BS_position;
        positions(cnt_tar+K,:) = tar_position;
        distance_tar = norm(tar_position-BS_position);
        distances(cnt_tar+K) = distance_tar;
        theta_tar = angle((tar_position-BS_position)*[1;1i]);
        theta_tars(cnt_tar) = theta_tar;
        thetas(cnt_tar+K) = theta_tar;
        beta_db_tar = generate_path_loss(distance_tar);
        betas(cnt_tar+K) = beta_db_tar;
    end
    
    sorted_theta_tars = sort(theta_tar);
    
    
    
elseif strcmp(UE_position_scenario,'random') && strcmp(target_position_scenario,'static')
    
    theta_tars = input_parameters.theta_tars;
    rho_tars = input_parameters.rho_tars;
    
        %% Generate UE positions    
    for cnt_ue = 1:K
        UE_position = (d_max_area/2-20).*rand(1,2); %x-y
        UE_position = UE_position.*(2.*randi([0,1],1,2)-1);
        UE_position = UE_position + BS_position;
        positions(cnt_ue,:) = UE_position;
        distance_UE = norm(UE_position-BS_position);
        distances(cnt_ue) = distance_UE;
        theta_UE = angle((UE_position-BS_position)*[1;1i]);
        thetas(cnt_ue) = theta_UE;
        beta_db_UE = generate_path_loss(distance_UE);
        betas(cnt_ue) = beta_db_UE;
    end



    %% Generate target positions

    thetas(K+1:end) = theta_tars;

    for cnt_tar = 1:T
        [tar_x, tar_y] = pol2cart(theta_tars(cnt_tar), rho_tars);
        target_position = [tar_x, tar_y] + [BS_position];
        positions(cnt_tar+K,:) = target_position;
        distance_tar = norm(target_position-BS_position);
        distances(cnt_tar+K) = distance_tar;
        beta_db_tar = generate_path_loss(distance_tar);
        betas(cnt_tar+K) = beta_db_tar;
    end    
       
elseif strcmp(UE_position_scenario,'static') && strcmp(target_position_scenario,'static')
    
    rho_UEs = input_parameters.rho_UEs;
    theta_UEs = input_parameters.theta_UEs;
    theta_tars = input_parameters.theta_tars;
    rho_tars = input_parameters.rho_tars;
    
        %% Generate UE positions    
    for cnt_ue = 1:K
        [UE_position_x, UE_position_y] = pol2cart(theta_UEs(cnt_ue), rho_UEs);
        UE_position = [UE_position_x, UE_position_y] + [BS_position];
       % UE_position = (d_max_area-10).*rand(1,2)+10; %x-y
        positions(cnt_ue,:) = UE_position;
        distance_UE = norm(UE_position-BS_position);
        distances(cnt_ue) = distance_UE;
        theta_UE = angle((UE_position-BS_position)*[1;1i]);
        thetas(cnt_ue) = theta_UE;
        beta_db_UE = generate_path_loss(distance_UE);
        betas(cnt_ue) = beta_db_UE;
    end



    %% Generate target positions

    thetas(K+1:end) = theta_tars;

    for cnt_tar = 1:T
        [tar_x, tar_y] = pol2cart(theta_tars(cnt_tar), rho_tars);
        target_position = [tar_x, tar_y] + [BS_position];
        positions(cnt_tar+K,:) = target_position;
        distance_tar = norm(target_position-BS_position);
        distances(cnt_tar+K) = distance_tar;
        beta_db_tar = generate_path_loss(distance_tar);
        betas(cnt_tar+K) = beta_db_tar;
    end

end


location_parameters.positions = positions;
location_parameters.distances = distances; 
location_parameters.thetas = thetas;
location_parameters.betas = betas;
end