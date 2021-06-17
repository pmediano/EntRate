%% var_to_tsdata
%
% Generate random multi-trial Gaussian VAR time series
%
% <matlab:open('var_to_tsdata.m') code>
%
%% Syntax
%
%     [X,E,mtrunc] = var_to_tsdata(A,V,m,N,mtrunc,decfac)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     V          residuals covariance matrix
%     m          number of observations per trial
%     N          number of trials (default: 1)
%     mtrunc     number of initial time observations to truncate or (default) empty for automatic calculation
%     decfac     initial transients decay factor (default: 1)
%
% _output_
%
%     X          multi-trial Gaussian VAR time series
%     E          residuals time series
%     mtrunc     actual number of initial time steps truncated
%
%% Description
%
% Return |N| time series of length |m| sampled from a VAR model with
% coefficients matrix |A|, and iid Gaussian residuals with covariance matrix
% |V|:
%
% <<eq_var.png>>
%
% If |mtrunc| is supplied it is taken to be the the number of initial
% (non-stationary transient) observations to truncate; otherwise (default) the
% spectral radius of |A| (see function <specnorm.html |specnorm|>) is
% calculated and used to estimate a suitable number |mtrunc| of observations to
% assumed stationarity (roughly till covariance decays to its stationary value).
% Set |decfac| > 1 for longer settle time.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
%% See also
%
% <specnorm.html |specnorm|> |
% <mvfilter.html |mvfilter|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [X,E,mtrunc] = var_to_tsdata(A,V,m,N,mtrunc,decfac)

if nargin < 4 || isempty(N), N = 1; end % single trial

if nargin < 5 || isempty(mtrunc) % automatic calclation - transients decay with rate given by VAR spectral radius
    if nargin < 6 || isempty(decfac), decfac = 1; end
    mtrunc = ceil(decfac*(-log(eps)+log(max(abs(eig(V)))))/(-log(specnorm(A))));
else
    assert(isscalar(mtrunc) && isint(mtrunc) && mtrunc >= 0,'truncation parameter must be a non-negative integer');
end

[VL,cholp] = chol(V,'lower');
assert(cholp == 0,'covariance matrix not positive-definite');

n = size(A,1);

if N > 1 % multi-trial

	E = zeros(n,m,N);
	for r = 1:N
		E(:,:,r) = VL*randn(n,mtrunc+m);
	end
	X = zeros(n,m,N);
	for r = 1:N
		X(:,:,r) = mvfilter([],A,E(:,:,r));
	end

else

	E = VL*randn(n,mtrunc+m);
	X = mvfilter([],A,E);

end

if mtrunc > 0
	X = X(:,mtrunc+1:mtrunc+m,:);
	if nargout > 1
		E = E(:,mtrunc+1:mtrunc+m,:);
	end
end
