%% mvgc_cdf
%
% Sample MVGC thoretical asymptotic cumulative distribution function
%
% <matlab:open('mvgc_cdf.m') code>
%
%% Syntax
%
%     P = mvgc_cdf(x,X,p,m,N,nx,ny,nz)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     x          vector of sample MVGC values
%     X          vector of actual MVGC values (zeros for null hypothesis of vanishing GC)
%     p          VAR model order
%     m          number of observations per trial
%     N          number of trials
%     nx         number of target ("to") variables
%     ny         number of source ("from") variables
%     nz         number of conditioning variables (default: 0)
%
% _output_
%
%     P          cumulative distribution probabilities evaluated at x
%
%% Description
%
% Return theoretical sample MVGC asymptotic cumulative distribution function for
% actual MVGCs in vector |X|, evaluated at values in vector |x|. For p-values,
% assume null hypothesis [H_0]: |X = 0| (see <mvgc_pval.html |mvgc_pval|>).
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] C. Ladroue, S. Guo, K. Kendrick and J. Feng, "Beyond element-wise
% interactions: Identifying complex interactions in biological processes",
% _PLoS ONE_ 4(9), 2009.
%
% [2] A. B. Barrett, L. Barnett and A. K. Seth, "Multivariate Granger causality
% and generalized variance", _Phys. Rev. E_ 81(4), 2010.
%
% [3] L. Barnett and A. K. Seth, "Behaviour of Granger causality under
% filtering: Theoretical invariance and practical application", _J. Neurosci.
% Methods_ 201(2), 2011.
%
%% See also
%
% <mvgc_cdfi.html |mvgc_cdfi|> |
% <mvgc_pval.html |mvgc_pval|> |
% <mvgc_confint.html |mvgc_confint|> |
% <mvgc_cval.html |mvgc_cval|> |
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function P = mvgc_cdf(x,X,p,m,N,nx,ny,nz,~)

xx = x; xx(isnan(x)) = [];
assert(all(xx >= 0), 'MVGC sample values must be non-negative');

XX = X; XX(isnan(X)) = [];
H0 = isempty(XX) || all(XX == 0);
if ~H0
	assert(all(XX >= 0),'MVGC actual values must be non-negative');
	if isscalar(X), X = X*ones(size(x)); else assert(isequal(size(X),size(x)),'MVGC actual values must match sample values'); end
end

if isempty(N), N = 1; end % single trial

if nargin < 8 || isempty(nz), nz = 0; end % unconditional

m = N*(m-p);                % effective number of observations (p-lag autoregression loses p observations per trial)
d = p*nx*ny;                % degrees of freedom = #{full model parameters} - #{reduced model parameters}
if H0
	P = chi2cdf(m*x,d);     % central
else
	P = ncx2cdf(m*x,d,m*X); % non-central approximation
end

