function numeric_correct = numeric_correct_complex(complex_number,sensitivity)
if abs(real(complex_number))<sensitivity
    real_part = 0;
else
    real_part = real(complex_number);
end

if abs(imag(complex_number))<sensitivity
    imag_part = 0;
else
    imag_part = imag(complex_number);
end

numeric_correct = real_part+1i*imag_part;

end