function H = Hurst_est(x,q)

% Lo, 1989

assert(iscolumn(x),'input must be a column vector');

n = length(x);

if nargin < 2 || isempty(q)
	q = 0;
else
	assert(q < n,'q too large');
end

x = x-sum(x)/n; % de-mean

y = cumsum(x);

Q = max(y)-min(y);

if q > 0 % Lo
	ssq = sum(x.*x);
	b = zeros(q,1);
	for j = 1:q
		i = (j+1):n;
		b(j) = sum(x(i).*x(i-j));
	end
	w = 1-(1:q)'/(q+1);
	ssq = ssq + 2*sum(w.*b);
	H = log(Q)/(sqrt(ssq/n)*log(n));
end
