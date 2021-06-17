function F = ss_to_pwcgc(A,C,K,V)

n = ss_parms(A,C,K,V);

F = nan(n);

KVL = K*chol(V,'lower');
KVK = KVL*KVL';
LV = log(diag(V));

for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

    [~,VR,rep] = ss2iss(A,C(r,:),KVK,V(r,r),K*V(:,r)); % "reduced" innovations covariance
    if sserror(rep,y), continue; end % check DARE report, bail out on error

    F(r,y) = log(diag(VR))-LV(r);
end
