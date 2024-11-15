function beta = generate_roundtrip_loss(carrier_freq, distance)
wavelength = 3e8/carrier_freq;
beta = wavelength^2/((4*pi)^3*(distance^4));
end