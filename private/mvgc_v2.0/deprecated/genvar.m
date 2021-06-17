function Y = genvar(A,X)

% Y = genvma(B,X)
%
% Generate vector autoregressive (VAR) process with coefficients A and input X.
%
% Implements:
%
%     Y(t) = X(t) + sum(k = 1:r) A(k)*Y(t-k)
%
% where r = min(p,t-1), with p the VAR model order.
%
% This function is essentially a wrapper for the MEX routine 'genvar_mex'.

global have_genvar_mex;

assert(isreal(X) && isa(X,'double') && ismatrix(X), 'Input must be a real, double-precision matrix')

[n,m] = size(X);

if isempty(A) % model order p = 0, just return input

    Y = X;

else

    assert(isreal(A) && isa(A,'double') && (ismatrix(A) || ndims(A) == 3), 'VAR coefficients must be a real, double-precision vector, matrix or 3D array')

	cvec = isvector(A); % vector of coefficients - apply filter to each row of X

	if ~cvec
		assert(size(A,1) == n && size(A,2) == n, 'VAR coefficient blocks must match input size');
	end

    if have_genvar_mex
        Y = genvar_mex(A,X,cvec);
    else
        Y = X;
        if cvec
			p = length(A);
			for t = 1:m
				for k = 1:p
					if k < t
						Y(:,t) = Y(:,t) + A(k)*Y(:,t-k);
					end
				end
			end
		else
			p = size(A,3);
			for t = 1:m
				for k = 1:p
					if k < t
						Y(:,t) = Y(:,t) + A(:,:,k)*Y(:,t-k);
					end
				end
			end
		end
    end

end
