function Y = genvma(B,X)

% Y = genvma(B,X,trunc)
%
% Generate vector moving-average (VMA) process with coefficients B and input X.
%
% Implements:
%
%     Y(t) = X(t) + sum(k = 1:r) B(k)*X(t-k)
%
% where r = min(q,t-1), with q the VMA model order.
%
% This function is essentially a wrapper for the MEX routine 'genvma_mex'.

global have_genvma_mex;

assert(isreal(X) && isa(X,'double') && ismatrix(X), 'Input must be a real, double-precision matrix')

[n,m] = size(X);

if isempty(B) % model order q = 0, just return input

	Y = X;

else

    assert(isreal(B) && isa(B,'double') && (ismatrix(B) || ndims(B) == 3), 'VMA coefficients must be a real, double-precision vector, matrix or 3D array')

	cvec = isvector(B); % vector of coefficients - apply filter to each row of X

	if ~cvec
		assert(size(B,1) == n && size(B,2) == n, 'VMA coefficient blocks must match input size');
	end

    if have_genvma_mex
        Y = genvma_mex(B,X,cvec);
    else
        Y = X;
        if cvec
			q = length(B);
			for t = 1:m
				for k = 1:q
					if k < t
						Y(:,t) = Y(:,t) + B(k)*X(:,t-k);
					end
				end
			end
		else
			q = size(B,3);
			for t = 1:m
				for k = 1:q
					if k < t
						Y(:,t) = Y(:,t) + B(:,:,k)*X(:,t-k);
					end
				end
			end
		end
    end

end
