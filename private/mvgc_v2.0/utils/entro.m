function [H,p] = entro(s,dim,nrm,tol)

assert(all(s(~isnan(s)) >= 0),'negative values found');

if nargin < 2 || isempty(dim)
	if isrow(s), dim = 2; else, dim = 1; end
end
if nargin < 3 || isempty(nrm), nrm = false;     end
if nargin < 4 || isempty(tol), tol = sqrt(eps); end

fac = nansum(s,dim);

p = bsxfun(@rdivide,s,fac);

H = -nansum(p.*log(p),dim);

H(fac < tol) = -Inf;

if nrm
	H = H./log(size(s,dim));
end
