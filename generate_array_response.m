function [array_derivative_response, array_response] = generate_array_response(array_geometry, number_of_antennas, DOA)
%%% array_geometry = 'ULA' (only option for now)
%%% DOA: [DOA_azimuth, DOA_elevation] 1x2 vector of DOA angles
if array_geometry == 'ULA'
    array_response = exp((1i)*pi*sin(DOA(1))*cos(DOA(2)).*(0:number_of_antennas-1)');
    array_der_coeff = ((1i)*pi*cos(DOA(1)).*(0:number_of_antennas-1)');
    array_derivative_response = array_der_coeff.*array_response;
end

%%%%%%%%%%%% Should there be (1/sqrt(number_of_antennas)) term???
end