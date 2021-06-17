function Y = mvfilter(B,A,X)

% Y = mvfilter(B,A,X)
%
% Implements a VARMA filter with VMA coefficients B, VAR coefficients A and input X.
%
% Implements:
%
%     Y(t) = X(t) + sum_k B(k)*X(t-k) + sum_k A(k)*Y(t-k)
%
% where the sums over k are for all indices for which the quantities are defined
%
% This function is a wrapper for the MEX routine 'mvfilter_mex'.

global have_mvfilter_mex;

assert(isreal(X) && isa(X,'double') && ismatrix(X), 'Input must be a real, double-precision matrix')

[n,m] = size(X);

if isempty(B)
	if isempty(A)
		Y = X;
		return
	end
	bvec = true;
	q = 0;
else
	assert(isreal(B) && isa(B,'double') && (ismatrix(B) || ndims(B) == 3), 'VMA coefficients must be a real, double-precision vector, matrix or 3D array')
	bvec = isvector(B); % vector of coefficients - apply filter to each row of X
	if bvec
		q = length(B);
	else
		assert(size(B,1) == n && size(B,2) == n, 'VMA coefficient blocks must match input size');
		q = size(B,3);
	end
end

if isempty(A)
	p = 0;
	avec = true;
else
	assert(isreal(A) && isa(A,'double') && (ismatrix(A) || ndims(A) == 3), 'VAR coefficients must be a real, double-precision vector, matrix or 3D array')
	avec = isvector(A); % vector of coefficients - apply filter to each row of X
	if avec
		p = length(A);
	else
		assert(size(A,1) == n && size(A,2) == n, 'VAR coefficient blocks must match input size');
		p = size(A,3);
	end
end

if have_mvfilter_mex
	Y = mvfilter_mex(B,bvec,q,A,avec,p,X);
	return
end

Y = X;
if bvec
	if avec % bvec, avec
		for t = 1:m
			for k = 1:q
				if k < t
					Y(:,t) = Y(:,t) + B(k)*X(:,t-k);
				end
			end
			for k = 1:p
				if k < t
					Y(:,t) = Y(:,t) + A(k)*Y(:,t-k);
				end
			end
		end
	else % bvec, ~avec
		for t = 1:m
			for k = 1:q
				if k < t
					Y(:,t) = Y(:,t) + B(k)*X(:,t-k);
				end
			end
			for k = 1:p
				if k < t
					Y(:,t) = Y(:,t) + A(:,:,k)*Y(:,t-k);
				end
			end
		end
	end
else % ~bvec
	if avec % ~bvec, avec
		for t = 1:m
			for k = 1:q
				if k < t
					Y(:,t) = Y(:,t) + B(:,:,k)*X(:,t-k);
				end
			end
			for k = 1:p
				if k < t
					Y(:,t) = Y(:,t) + A(k)*Y(:,t-k);
				end
			end
		end
	else % ~bvec, ~avec
		for t = 1:m
			for k = 1:q
				if k < t
					Y(:,t) = Y(:,t) + B(:,:,k)*X(:,t-k);
				end
			end
			for k = 1:p
				if k < t
					Y(:,t) = Y(:,t) + A(:,:,k)*Y(:,t-k);
				end
			end
		end
	end
end
