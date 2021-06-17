function F = ss_to_cwiogc(A,C,K,V,inout)

% in/out group GCs

[n,r] = ss_parms(A,C,K,V);

switch lower(inout)
	case 'in' , gcin = true;
	case 'out', gcin = false;
	otherwise , error('in/out parameter must be ''in'' or ''out''');
end

F = nan(n,1);

KVL = K*chol(V,'lower');
KVK = KVL*KVL';

for i = 1:n
	if gcin
		r = i;
	else
		r = 1:n; r(i) = [];
	end

    [~,VR,rep] = ss2iss(A,C(r,:),KVK,V(r,r),K*V(:,r)); % reduced model innovations covariance
    if sserror(rep,i), continue; end % check DARE report, bail out on error

	F(i) = logdet(VR) - logdet(V(r,r));
end
