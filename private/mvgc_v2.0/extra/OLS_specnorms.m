function rho = OLS_specnorms(X,p,verb)

[n,m,N] = size(X);

assert(isscalar(p) && isint(p) && p > 0 && p < m,'maximum model order must be a positive integer less than the number of observations');

if nargin < 3 || isempty(verb), verb = false; end

p1 = p+1;
pn = p*n;
p1n = p1*n;

I = eye(n);

X = demean(X); % no constant term, normalize

% store lags

XX = zeros(n,p1,m+p,N);
for k = 0:p
    XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
end

rho = nan(p,1); % spectral norm

% loop through model orders

for k = 1:p

	if verb, fprintf('model order = %2d : ',k); end

	k1 = k+1;
	kn1 = (k-1)*n;
	M  = N*(m-k);

	X0 = reshape(XX(:,1,   k1:m,:),n,  M);
	XL = reshape(XX(:,2:k1,k1:m,:),n*k,M);

	wstate = warning('off','all');
	lastwarn('');
	A = X0/XL;                  % OLS (QR decomposition)
	wmsg = lastwarn;
	warning(wstate);
	if ~isempty(wmsg) % rank-deficient? Not necessarily a show-stopper
		if ~verb, fprintf('model order = %d',k); end
		fprintf(2,'WARNING: VAR estimation may be problematic (%s)',wmsg);
		% not necessarily a show-stopper - carry on
	end
	if isbad(A)       % something went badly wrong; show stopper
		if ~verb, fprintf('model order = %d',k); end
		fprintf(2,'WARNING: VAR estimation failed\n');
		continue % show-stopper
	end

	% calculate spectral radius

	AA = [reshape(A,n,k*n); eye(kn1) zeros(kn1,n)];
	rho(k) = max(abs(eig(AA,'nobalance')));


	if verb, fprintf('rho = %6.4f\n',rho(k)); end
end
