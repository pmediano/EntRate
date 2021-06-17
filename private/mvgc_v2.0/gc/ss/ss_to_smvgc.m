function f = ss_to_smvgc(A,C,K,V,x,y,fres)

[n,r] = ss_parms(A,C,K,V);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
w = [y z];
xr = 1:length(x);

h = fres+1;
f = nan(1,h);

H   = ss2trfun(A,C,K,fres);
VL  = chol(V,'lower');
PVL = chol(parcov(V,w,x),'lower');

% Note: in theory we shouldn't have to take the real part of the determinants of
% the (Hermitian, positive-definite) matrices in the calculation of the f(k),
% since they should be real. However, Matlab's det() will sometimes return a
% tiny imaginary part.

if isempty(z) % unconditional (note: does not require reduced model)

    for k = 1:h
        HVL  = H(x,:,k)*VL;
        SR   = HVL*HVL';
        HR   = H(x,y,k)*PVL;
        f(k) = logdet(SR) - logdet(SR-HR*HR');
    end

else % conditional

    KVL = K*VL;
    CR  = C(r,:);
    [KR,SIGR,rep]  = ss2iss(A,CR,KVL*KVL',V(r,r),K*V(:,r)); % reduced model Kalman gain and innovations covariance
    if sserror(rep), return; end % check DARE report, bail out on error

    BR    = ss2itrfun(A,CR,KR,fres);
    SR    = SIGR(xr,xr); % reduced model spectrum is flat!
    LDSR  = logdet(SR);
    for k = 1:h
        HR   = BR(xr,:,k)*H(r,w,k)*PVL;
        f(k) = LDSR - logdet(SR-HR*HR');
    end

end
