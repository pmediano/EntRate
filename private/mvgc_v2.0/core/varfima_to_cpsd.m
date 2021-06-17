function [S,H] = varfima_to_cpsd(A,B,fiparms,V,fres)

have_var  = ~isempty(A);
have_vma  = ~isempty(B);
fint      = ~isempty(fiparms);

if fint
	assert(isvector(fiparms) && length(fiparms) == 3,'fractional integration parameters must be a 2-vector [d,r]');
	d = fiparms(1);
	r = fiparms(2);
	arform = fiparms(3);
end

assert(ismatrix(V),'covariance must be a square matrix');
[n,n1] = size(V);
assert(n1 == n,'covariance matrix not square');

S = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')
H = NaN; % ensure a "bad" return value if anything goes wrong (see routine 'isbad')

[L,cholp] = chol(V,'lower');
if cholp, return; end % show stopper

I = eye(n);

h = fres+1; % all calculations over [0,2*pi)

if have_var
	[n1,n2,~] = size(A);
	assert(n1 == n && n2 == n,'VAR coefficients array does not match covariance matrix');
	AF = fft(cat(3,I,-A),2*fres,3);
end

if have_vma
	[n1,n2,~] = size(B);
	assert(n1 == n && n2 == n,'VMA coefficients array does not match covariance matrix');
	BF = fft(cat(3,I,B),2*fres,3);
end

if fint
	c = fracint_coeffs(d,r,arform);
	if arform
		bf = 1./fft(c,2*fres);
	else
		bf = fft(c,2*fres);
	end
end

H = zeros(n,n,h);
if have_var
	if have_vma
		if fint
			for k = 1:h
				H(:,:,k) = bf(k)*(AF(:,:,k)\BF(:,:,k));
			end
		else
			for k = 1:h
				H(:,:,k) = AF(:,:,k)\BF(:,:,k);
			end
		end
	else
		if fint
			for k = 1:h
				H(:,:,k) = bf(k)*(AF(:,:,k)\I);
			end
		else
			for k = 1:h
				H(:,:,k) = AF(:,:,k)\I;
			end
		end
	end
else
	if have_vma
		if fint
			for k = 1:h
				H(:,:,k) = bf(k)*BF(:,:,k);
			end
		else
			H = BF;
		end
	else
		if fint
			for k = 1:h
				H(:,:,k) = bf(k);
			end
		else
			for k = 1:h
				H(:,:,k) = I;
			end
		end
	end
end

S = zeros(n,n,h);
for k = 1:h
	HLk = H(:,:,k)*L;
	S(:,:,k) = HLk*HLk';
end
