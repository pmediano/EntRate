function c = fracint_coeffs(d,r,arform)

% Fractional integration coefficients in AR or MA form

if arform % return AR coefficients
	dd = 1+d;
else      % return MA coefficients
	dd = 1-d;
end

c = zeros(1,r+1);
c(1) = 1;
for k = 1:r
	c(k+1) = (1-dd/k)*c(k);
end
