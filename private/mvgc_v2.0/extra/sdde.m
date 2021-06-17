function y = sdde(a,c,x,dt,use_sdde_mex)

% Naive linear Stochastic Delay-Differential Equation simulation

[n,m] = size(x);

assert(isvector(a) && length(a) == n);

[C,r] = size(c);

assert(ismatrix(c) && r == 4);

% Scale Wiener noise

x = sqrt(dt)*x;

% Setup/scale decay (AR) coefficients

a = 1-dt./a(:);

% Setup/scale delay-connection parameters

dt

i = c(:,1);
j = c(:,2);
s = dt./c(:,3);
d = round(c(:,4)/dt);
assert(all(isint(i) & i >= 1 & i <= n));
assert(all(isint(j) & j >= 1 & j <= n));
assert(all(i ~= j));
assert(all(d > 0));
a
d
s

if use_sdde_mex
	y = sdde_mex(a,uint64(i-1),uint64(j-1),s,uint64(d),x);
else
	y = x;
	for t = 2:m
		y(:,t) = y(:,t) + a.*y(:,t-1);
		for k = 1:C
			if t > d(k)
				y(i(k),t) = y(i(k),t) + s(k)*y(j(k),t-d(k));
			end
		end
	end
end
