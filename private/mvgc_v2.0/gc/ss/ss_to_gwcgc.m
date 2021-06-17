function F = ss_to_gwcgc(A,C,K,V,group)

% inter-group (conditional) GCs

[n,r] = ss_parms(A,C,K,V);

assert(iscell(group) && isrow(group),     'groups must be specified as a cell (row) vector of indices');
idx = cell2mat(group);
assert(length(unique(idx)) == length(idx),'groups indices must be unique and non-overlapping');
assert(all(idx >=1 & idx <= n),           'some group indices out of range');

g = length(group);

F = NaN(g);

KVL = K*chol(V,'lower');
KVK = KVL*KVL';

for b = 1:g
	y = group{b};
    r = 1:n; r(y) = []; % omit group b

    [~,VR,rep] = ss2iss(A,C(r,:),KVK,V(r,r),K*V(:,r)); % "reduced" innovations covariance
    if sserror(rep,b), continue; end % check DARE report, bail out on error

    for a = 1:g
        if a == b, continue; end
        x = group{a};
        xr = findin(x,r); % indices of group{a} in r
        F(a,b) = logdet(VR(xr,xr))-logdet(V(x,x));
    end
end
