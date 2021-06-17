function [Y,nouts] = routl(X,sdfac,orep)

% Replace outliers in each variable in matrix X (variables x observations x
% trials). Outliers are defined as less than (mean - sdfac*sd) or greater than
% (mean + sdfac*sd) where sd is the standard deviation of X.
%
% Replacement strategy depends on 'orep':
%
% 1 : replace with mean of non-outliers
% 2 : replace with random non-outlier
% 3 : replace with Gaussian deviate with mean & std. dev of non-outliers

if sdfac == 0 % do nothing
    Y = X;
    nouts = 0;
    return;
end

assert(orep == 1 || orep == 2 || orep == 3,'''orep'' must be 1, 2, or 3');

[n,m,N] = size(X);
X = X(:,:);

Y = X;
nouts = zeros(n,1);

for i = 1:n
    xmean      = mean(X(i,:));
    xsdev      = std(X(i,:));
    outs       = X(i,:) < xmean - sdfac*xsdev | X(i,:) > xmean + sdfac*xsdev; % logical idx of outliers
    nouts(i)   = nnz(outs);  % number of outliers
	Z          = X(i,~outs); % non-outliers
	if     orep == 1 % replace with mean of non-outliers
        Y(i,outs)  = mean(Z)*ones(1,nouts(i));
    elseif orep == 2 % replace with random non-outlier
        Y(i,outs)  = Z(randi(length(Z),[1 nouts(i)]));
	elseif orep == 3 % replace with Gaussian with mean & std. dev of non-outliers
        Y(i,outs)  = mean(Z) + std(Z)*randn(1,nouts(i));
	end
end

if N > 1 % multi-trial
    Y = reshape(Y,[n m N]);
end
