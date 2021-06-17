function Y = gdfilter(X,parms)

% Generalised differencing of matrix. Calculates:
%
%     Y(t) = (1-L)^(-d) X(t)
%
% where L is the lag (backshift) operator. The differencing order
% d may be fractional. Differencing is applied row-wise to X.
%
% INPUTS
%
% X        input matrix
% parms    parameter 3-vector or 3-cell vector containing:
%          d      - differencing order (can be fractional)
%          r      - number of AR/MA coefficients to calculate
%          arform - use AR, rather than MA filter
%
% OUTPUTS
%
% Y        output matrix

assert(isvector(parms) && length(parms) == 3,'parameters must be a 3-vector or 3-cell vector');
if iscell(parms)
	d = parms{1};
	r = parms{2};
	arform = parms{3};
	dvec = ~isscalar(d);
else
	d = parms(1);
	r = parms(2);
	arform = logical(parms(3));
	dvec = false;
end

if dvec
	[n,m] = size(X);
	assert(length(d) == n,'differencing exponent vector doesn''t match input vector');
	if isscalar(r)
		r = r*ones(n,1);
	else
		assert(length(r) == n,'number of coefficients vector doesn''t match input vector');
	end
	if isscalar(arform)
		arform = arform&true(n,1);
	else
		assert(length(arform) == n,'AR form vector doesn''t match input vector');
	end
	Y = zeros(n,m);
	for i = 1:n
		if abs(d(i)) > eps
			c = fracint_coeffs(d(i),r(i),arform(i));
			if arform % use AR form
				Y(i,:) = mvfilter([],-c(2:end),X(i,:));
			else      % use MA form
				Y(i,:) = mvfilter(c(2:end),[],X(:,i));
			end
		else
			Y = X;
		end
	end
else
	if abs(d) > eps
		c = fracint_coeffs(d,r,arform);
		if arform % use AR form
			Y = mvfilter([],-c(2:end),X);
		else      % use MA form
			Y = mvfilter(c(2:end),[],X);
		end
	else
		Y = X;
	end
end
