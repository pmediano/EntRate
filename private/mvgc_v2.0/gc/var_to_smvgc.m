function f = var_to_smvgc(A,V,x,y,fres)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
w = [y z];

h = fres+1;
f = nan(1,h);

H   = var2trfun(A,fres);
VL  = chol(V,'lower');
PVL = chol(parcov(V,w,x),'lower');

if isempty(z) % unconditional (note: does not require reduced model)

    for k = 1:h
        HVL  = H(x,:,k)*VL;
        SR   = HVL*HVL';
        HR   = H(x,y,k)*PVL;
        f(k) = logdet(SR) - logdet(SR-HR*HR');
    end

else % conditional

	nx = length(x);
	ny = length(y);
	nz = length(z);
	nr = nx+nz;
	xr = 1:nx;

	pn   = p*n;
	pn1  = pn-n;
	pny  = p*ny;
	pny1 = pny-ny;

	% Solve the shrunken DARE

	[KT,VR,rep] = var2riss(A,V,y,r);
    if sserror(rep), return; end % check DARE report, bail out on error

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

	% Calculate spectral GC

    SR    = VR(xr,xr); % reduced model spectrum is flat!
    LDSR  = logdet(SR);
    for k = 1:h
        HR   = BR(xr,:,k)*H(r,w,k)*PVL;
        f(k) = LDSR - logdet(SR-HR*HR');
    end

end
