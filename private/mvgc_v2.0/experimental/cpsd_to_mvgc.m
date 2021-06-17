function F = cpsd_to_mvgc(S,V,x,y,maxiters,numtol,verb)

if nargin < 5, maxiters = []; end % cpsd_to_var default
if nargin < 6, numtol   = []; end % cpsd_to_var default
if nargin < 7, verb     = []; end % cpsd_to_var default

n = size(S,1);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];
xr = 1:length(x);        % index of x in reduced quantities

F = NaN;

[~,VR,converged,~,mre] = cpsd_to_var(S(r,r,:),maxiters,numtol,verb);
if ~converged
    fprintf(2,'WARNING: spectral factorisation failed to converge in %d iterations (relative residual = %e)\n',mre);
end

F = logdet(VR(xr,xr)) - logdet(V(x,x));
