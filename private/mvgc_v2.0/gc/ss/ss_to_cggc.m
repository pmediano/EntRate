function F = ss_to_cggc(A,C,K,V,x)

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

F = 0;

KVL = K*chol(V,'lower');
KVK = KVL*KVL';

for i = 1:nx
    r = [x(i) z];

    [~,VR,rep] = ss2iss(A,C(r,:),KVK,V(r,r),K*V(:,r)); % reduced model innovations covariance
    if sserror(rep,i), continue; end % check DARE report, bail out on error

    F = F + log(VR(1,1)) - log(V(x(i),x(i)));
end
