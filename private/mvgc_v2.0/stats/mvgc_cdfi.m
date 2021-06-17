%% mvgc_cdfi
%
% Sample MVGC thoretical asymptotic inverse cumulative distribution function
%
% <matlab:open('mvgc_cdfi.m') code>
%
%% Syntax
%
%     x = mvgc_cdfi(P,X,p,m,N,nx,ny,nz)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     P          vector of MVGC cumulative distribution probabilities
%     X          vector of actual MVGC values
%     p          VAR model order
%     m          number of observations per trial
%     N          number of trials
%     nx         number of target ("to") variables
%     ny         number of source ("from") variables
%     nz         number of conditioning variables (default: 0)
%
% _output_
%
%     x          vector of MVGC values
%
%% Description
%
% Return theoretical sample MVGC asymptotic inverse cumulative distribution
% function for actual MVGCs in vector |X|, evaluated at probabilities in vector
% P. To calculate the critical MVGC value for a significance level |alpha|,
% assume null hypothesis [H_0]: |X = 0| and set |P = 1-alpha| (see
% <mvgc_cval.html |mvgc_cval|>). For confidence intervals at level |alpha|, set
% |X| to the estimated MVGC and set |P = alpha| for lower bound and |P =
% 1-alpha| for upper bound (see <mvgc_confint.html |mvgc_confint|>).
%
% For more details see <mvgc_cdf.html |mvgc_cdf|>.
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
% <mvgc_cdf.html |mvgc_cdf|> |
% <mvgc_pval.html |mvgc_pval|> |
% <mvgc_confint.html |mvgc_confint|> |
% <mvgc_cval.html |mvgc_cval|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function x = mvgc_cdfi(P,X,p,m,N,nx,ny,nz,~)

PP = P; PP(isnan(P)) = [];
assert(all(PP >= 0) && all(PP <= 1), 'probabilities must lie between 0 and 1');

XX = X; XX(isnan(X)) = [];
H0 = isempty(XX) || all(XX == 0);
if ~H0
	assert(all(XX >= 0),'MVGC actual values must be non-negative');
	if isscalar(X), X = X*ones(size(P)); else assert(isequal(size(X),size(P)),'MVGC actual values must match probabilities'); end
end

if isempty(N), N = 1; end % single trial

if nargin < 8 || isempty(nz), nz = 0; end % unconditional

m = N*(m-p);                % effective number of observations (p-lag autoregression loses p observations per trial)
d = p*nx*ny;                % degrees of freedom = #{full model parameters} - #{reduced model parameters}
if H0
	x = chi2inv(P,d)/m;     % central
else
	x = ncx2inv(P,d,m*X)/m; % non-central approximation
end
