function [A_matrix, A_derivative_matrix] = generate_A_matrices(number_of_antennas, AOD)
%%% array_geometry = 'ULA' (only option for now)
%%% AOD: [AOD_azimuth, AOD_elevation] 1x2 vector of DOA angles

x = exp((1i)*pi*sin(AOD(1))*cos(AOD(2)));
der_term = (1i)*cos(AOD(1))*pi;
for m_1= 1:number_of_antennas
    for m_2 = 1:number_of_antennas
        A_matrix(m_1,m_2) = x^(m_1+m_2-2);
        A_derivative_matrix(m_1,m_2) = (m_1+m_2-2)*der_term*(x^(m_1+m_2-2));
    end
end


end
