function e = nanequal(x,y,tol)

exact = nargin < 3 || isempty(tol);

e = false;
if ~isequal(size(x),size(y)), return; end

inx = isnan(x);
iny = isnan(y);
if ~isequal(inx,iny), return; end

x(inx) = [];
y(iny) = [];
if exact
	if ~isequal(x,y), return; end
else
	if any(abs(x(:)-y(:)) > tol), return; end
end

e = true;
