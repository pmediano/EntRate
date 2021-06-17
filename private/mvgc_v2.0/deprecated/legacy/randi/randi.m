function r = randi(imax,n,m)

% Amended v2.0 to take 3 arguments

%fprintf(2,'legacy imax!\n');
assert(nargin > 0,'too few arguments');
if     nargin == 1
    r = 1+floor(imax*rand);
elseif nargin == 2
    r = 1+floor(imax*rand(n));
elseif nargin == 3
    r = 1+floor(imax*rand(n,m));
end
