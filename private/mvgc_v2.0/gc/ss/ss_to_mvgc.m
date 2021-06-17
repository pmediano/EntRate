%IMPORTANT: You MUST run ss_info before using this function!

function F = ss_to_mvgc(A,C,K,V,x,y)

n = ss_parms(A,C,K,V);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
xr = 1:length(x);        % index of x in reduced quantities

F = NaN;

KVL = K*chol(V,'lower');

[~,VR,rep] = ss2iss(A,C(r,:),KVL*KVL',V(r,r),K*V(:,r)); % reduced model innovations covariance
if sserror(rep), return; end % check DARE report, bail out on error

F = logdet(VR(xr,xr)) - logdet(V(x,x));
