function f = var_to_sgwggc(A,V,group,fres)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

g = check_group(group,n);

pn   = p*n;
pn1  = pn-n;

h = fres+1;
f = nan(g,h);

H = var2trfun(A,fres);

for a = 1:g
	x = group{a};
    z = 1:n; z(x) = []; % omit x
    nx = length(x);
	ny = nx-1;
	nz = length(z);
	nr = nz+1;
fprintf('\n\t\tG : group %d of %d : %2d sources ',a,g,nx);

	pny  = p*ny;
	pny1 = pny-ny;

	fa = zeros(1,h);
	sserr = false;
	for i = 1:nx
fprintf('.');
		xi = x(i);
		y = x; y(i) = [];
		r = [xi z];

		[KT,VR,rep] = var2riss(A,V,y,r);
		if sserror(rep,i) % check DARE report, bail out on error
			sserr = true;
			break;
		end

		% Calculate reduced SS parameters from shrunken DARE (note: VR is the same)

		AR = [reshape(A,n,pn); eye(pn1) zeros(pn1,n)];
		CR = reshape(A(r,:,:),nr,pn);

		%KR = kt2kr_mex(KT,r,y); % mex version of below: actually, no faster!

		KR = zeros(pn,nr);
		KR(r,:) = eye(nr);
		kn = 0;
		for ky = 0:ny:pny1
			KR(kn+y,:) = KT(ky+1:ky+ny,:);
			kn = kn+n;
		end
		BR = ss2itrfun(AR,CR,KR,fres);

		SR = VR(1,1); % reduced model spectrum is flat!
		LSR = log(SR);

		w = 1:n; w(xi) = [];
		PVL = chol(parcov(V,w,xi),'lower');

		for k = 1:h
			HR   = BR(1,:,k)*H(r,w,k)*PVL;
			fa(k) = fa(k) + LSR - log(SR-HR*HR');
		end
	end
	if sserr, continue; end % bail out on earlier error
	f(a,:) = fa;
end
