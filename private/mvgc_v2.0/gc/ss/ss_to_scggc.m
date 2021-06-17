function f = ss_to_scggc(A,C,K,V,x,fres)

[n,r] = ss_parms(A,C,K,V);

if nargin < 5 || isempty(x)
    x = 1:n;
    nx = n;
    z = [];
else
    x = x(:)'; % vectorise
	assert(length(unique(x)) == length(x),'x indices must be unique');
    assert(all(x >=1 & x <= n),'some x indices out of range');
    nx = length(x);
    z = 1:n; z(x) = []; % omit x
end

h = fres+1;
f = zeros(h,1);

H   = ss2trfun(A,C,K,fres);
KVL = K*chol(V,'lower');
KVK = KVL*KVL';

for i = 1:nx
    r = [x(i) z];
    CR = C(r,:);
    [KR,VR,rep] = ss2iss(A,CR,KVK,V(r,r),K*V(:,r)); % reduced model innovations covariance
    if sserror(rep,i), continue; end % check DARE report, bail out on error

    BR = ss2itrfun(A,CR,KR,fres);
    SR = VR(1,1); % reduced model spectrum is flat!
    LSR = log(SR);

    w = 1:n; w(x(i)) = [];
    PVL = chol(parcov(V,w,x(i)),'lower');

    for k = 1:h
        HR   = BR(1,:,k)*H(r,w,k)*PVL;
        f(k) = f(k) + LSR - log(SR-HR*HR');
    end
end
