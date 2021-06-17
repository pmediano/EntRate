function f = ss_to_sgwcgc(A,C,K,V,group,fres)

% inter-group (conditional) spectral GCs

[n,r] = ss_parms(A,C,K,V);

assert(iscell(group) && isrow(group),     'groups must be specified as a cell (row) vector of indices');
idx = cell2mat(group);
assert(length(unique(idx)) == length(idx),'groups indices must be unique and non-overlapping');
assert(all(idx >=1 & idx <= n),           'some group indices out of range');

g = length(group);

h = fres+1;
f = nan(g,g,h);

H   = ss2trfun(A,C,K,fres);
KVL = K*chol(V,'lower');
KVK = KVL*KVL';

% for efficiency (in time if not memory) we pre-compute partial covariances

PVL = cell(g,1);
for a = 1:g
    x = group{a};
    w = 1:n; w(x) = []; % omit group a
    PVL{a} = chol(parcov(V,w,x),'lower'); % pre-compute the partial covariances for efficiency
end

for b = 1:g
    r = 1:n; r(group{b}) = []; % omit group b

    CR = C(r,:);
    [KR,VR,rep] = ss2iss(A,CR,KVK,V(r,r),K*V(:,r)); % "reduced" innovations covariance
    if sserror(rep,b), continue; end % check DARE report, bail out on error

    BR = ss2itrfun(A,CR,KR,fres);

    for a = 1:g
        if a == b, continue; end
        x = group{a};
        xr = findin(x,r);   % indices of group{a} in r
        w = 1:n; w(x) = []; % omit group a

        SR  = VR(xr,xr);   % reduced model spectrum is flat!
        LDSR = logdet(SR);

        PVLa = PVL{a};
        for k = 1:h
            HR = BR(xr,:,k)*H(r,w,k)*PVLa;
            f(a,b,k) = LDSR - logdet(SR-HR*HR');
        end
    end
end
