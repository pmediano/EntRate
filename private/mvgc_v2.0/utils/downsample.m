function [Y,fs] = downsample(X,k,FS,o)

% Downsample time series data in X. The downsample factor k must be a positive
% integer less than or equal to the number of observations. The offset o must be
% a non-negative integer less than k, or -1 (default is o = 0). If o == -1, then
% for offset zero to k-1, the corresponding downsampled series are stacked (the
% result will always be multi-trial).

if nargout > 1
	assert(nargin > 2 && ~isempty(FS),'Sampling rate must be specified if new sampling rate requested');
	fs = FS/k;
end

if nargin < 4 || isempty(o), o = 0; end

[n,m,N] = size(X);

assert(isscalar(k) && isint(k) && k >= 0 && k <= m,              'Bad downsample factor (must be >= 1 and <= number of observations)');
assert(isscalar(o) && isint(o) && ((o >= 0 && o < k) || o == -1),'Bad offset (must be non-negative and < downsample factor, or equal to -1)');

if k == 1; % nothing to do!
    Y = X;
    return
end

if o == -1
    mk = floor(m/k);
    m = k*mk;
    Y = zeros(n,mk,N*k);
    for o = 0:k-1 % for each offset
        Y(:,:,N*o+1:N*(o+1)) = X(:,o+1:k:m,:);
    end
else
    Y = X(:,o+1:k:m,:);
end
