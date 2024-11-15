function beta_db = generate_path_loss(distance)
%%% Generates path loss terms for sub6GHz band according to the 3GPP Urban
%%% Microcell Model
beta_db = -30.5 -36.7*log10(distance);
end