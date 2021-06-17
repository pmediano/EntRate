function f = var_to_scwiogc(A,V,inout,fres)

% in/out group spectral GCs

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

switch lower(inout)
	case 'in',  gcin = true;
	case 'out', gcin = false;
	otherwise, error('in/out parameter must be ''in'' or ''out''');
end

h = fres+1;
f = nan(n,h);

H   = var2trfun(A,fres);
VL  = chol(V,'lower');

if gcin
	for i = 1:n
		x = i;
		y = 1:n; y(x) = [];

		PVL = chol(parcov(V,y,x),'lower');

		for k = 1:h
			HVL  = H(x,:,k)*VL;
			SR   = HVL*HVL';
			HR   = H(x,y,k)*PVL;
			f(i,k) = log(SR) - log(SR-HR*HR'); % 1-dim (no determinant required!)
		end
	end
else
	for i = 1:n
		y = i;
		x = 1:n; x(y) = [];

		PVL = chol(parcov(V,y,x),'lower');

		for k = 1:h
			HVL  = H(x,:,k)*VL;
			SR   = HVL*HVL';
			HR   = H(x,y,k)*PVL;
			f(i,k) = logdet(SR) - logdet(SR-HR*HR');
		end
	end
end
