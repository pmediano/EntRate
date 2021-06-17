function f = var_to_sgwiogc(A,V,group,inout,fres)

% in/out group spectral GCs

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(group,n);

switch lower(inout)
	case 'in',  gcin = true;
	case 'out', gcin = false;
	otherwise, error('in/out parameter must be ''in'' or ''out''');
end

h = fres+1;
f = nan(g,h);

H   = var2trfun(A,fres);
VL  = chol(V,'lower');

for a = 1:g
	if gcin
		x = group{a};
		y = 1:n; y(x) = [];
	else
		y = group{a};
		x = 1:n; x(y) = [];
	end

	PVL = chol(parcov(V,y,x),'lower');

    for k = 1:h
        HVL  = H(x,:,k)*VL;
        SR   = HVL*HVL';
        HR   = H(x,y,k)*PVL;
        f(a,k) = logdet(SR) - logdet(SR-HR*HR');
    end
end
