%% cpsd_to_var
%
% Spectral factorisation: calculate VAR parameters from cross-power spectral density
%
% <matlab:open('cpsd_to_var.m') code>
%
%% Syntax
%
%     [H,V,iters] = cpsd_to_var(S,G0,maxiters,numtol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     S          cross-power spectral density (cpsd) matrix
%     maxiters   maximum iterations (default: 100)
%     numtol     numerical tolerance (default: 1e-9)
%     verb       verbosity flag (default: 0)
%
% _output_
%
%     H          VAR transfer function matrix
%     V          VAR residuals covartiance matrix
%     converged  Boolean flag to indicate whether the process converged to the specified tolerance
%     iters      number of iterations performed
%     mre        maximum relatabsolute elemnt-wise error of factored CPSD with respect to supplied CPSD
%
%% Description
%
% Calculates the transfer function |H| and residuals covariance matrix |V|
% from the cpsd |S| and covariance matrix |G0| of a VAR process, using Wilson's
% spectral factorisation algorithm [2]. If |G0| is not supplied (default), then
% |G| is calculated from |S| (see <cpsd_to_autocov.html |cpsd_to_autocov|>) and
% |G0| set to |G(:,:,1)|. The actual number of iterations performed is returned
% in |iters|. If the algorithm fails to converge to numerical tolerance |numtol|
% within |maxiters| iterations, an exception |MVGC:XMaxItrs| is thrown.
%
% *_Note:_* to calculate the VAR coefficients, use the utility function
% <trfun2var.html |trfun2var|>.
%
% Adapted from original code with the kind permission of
% <http://math.iisc.ernet.in/~rangaraj/ Prof. G. Rangarajan> of the Dept. of
% Mathematics, Indian Institute of Science, Bangalore, India; see [3,4] for
% applications to Granger-causal analysis.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] G. T. Wilson, "The factorization of matricial spectral densities", _SIAM
% Journal on Applied Mathematics_, 23(4), 1972.
%
% [3] M. Dhamala, G. Rangarajan and M. Ding, "Estimating Granger causality from
% Fourier and wavelet transforms of time series data", _Phys. Rev. Lett._ 100,
% 2008.
%
% [4] M. Dhamala, G. Rangarajan and M. Ding, "Analyzing information flow in
% brain networks with nonparametric Granger causality", _NeuroImage_ 41, 2008.
%
%% See also
%
% <cpsd_to_autocov.html |cpsd_to_autocov|> |
% <autocov_to_var.html |autocov_to_var|> |
% <trfun2var.html |trfun2var|>
%
%% Copyright notice
%
% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file
% <matlab:open('license.txt') license.txt> in root directory for licensing
% terms._
%
%%

function [H,V,converged,iters,mre] = cpsd_to_var(S,maxiters,numtol,verb)

if nargin < 3 || isempty(numtol),   numtol   = 1e-9; end
if nargin < 2 || isempty(maxiters), maxiters = min(500,floor(sqrt(10/numtol))); end
if nargin < 4 || isempty(verb),     verb     = 0;    end

[n,~,h] = size(S);
h2 = 2*(h-1);

SX = cat(3,S,conj(S(:,:,h-1:-1:2))); % extend CPSD

% initialise P0
gamma =  ifft(SX,[],3);
gamma0 = gamma(:,:,1);
gamma0 = real((gamma0+gamma0')/2); % enforce symmetry, zero out imaginaries on diagonal
P0 = chol(gamma0);

P = repmat(P0,[1,1,h]);
P = cat(3,P,conj(P(:,:,h-1:-1:2)));

I = eye(n);

U = zeros(size(SX));
for i = 1:h2
    U(:,:,i) = chol(SX(:,:,i),'lower');
end

g = zeros(n,n,h2);
SF = zeros(n,n,h);

for iters = 1:maxiters

    if verb, fprintf('iteration %2d ...',iters); end

    for i = 1:h2
        % Equivalent to: g(:,:,i) = P(:,:,i)\SX(:,:,i)/P(:,:,i)' + I;
        C = P(:,:,i)\U(:,:,i);
        g(:,:,i) = C*C' + I;
    end

    % []+ operation
    gam = real(ifft(g,[],3));
	gam(:,:,1) = gam(:,:,1)/2; % take half of the zero lag
	gp0 = gam(:,:,1);
	gam(:,:,h+1:end) = 0;      % zero out negative powers.
	gp = fft(gam,[],3);        % reconstitute

    T = -tril(gp0,-1);
    T = T-T';

    P_prev = P;
    for i = 1:h2
        P(:,:,i) = P(:,:,i)*(gp(:,:,i) + T);
    end

    P0_prev = P0;
    P0 = P0*(gp0 + T);

    % Relative Cauchy error. Check on S is expensive, so check P0 first, then P and only then S
    converged = false;
	mre = maxrelerr(P0_prev,P0);
    if mre < numtol
		mre = maxrelerr(P_prev,P);
        if mre < numtol
            for i = 1:h
                SF(:,:,i) = P(:,:,i)*P(:,:,i)';
            end
            mre = maxrelerr(S,SF);
            converged = mre < numtol;
        end
    end

    if verb, fprintf(' mre = %.2e\n',mre); end

    if converged
        break
    end
end

H = zeros(n,n,h);
for i = 1:h
    H(:,:,i) = P(:,:,i)/P0;
end

V = P0*P0';

%% Maximum absolute element-wise error of y with respect to x

function mre = maxrelerr(x,y)

d = abs(y-x);
x = abs(x);
x(x <= 2*eps) = 1;
e = d./x;
mre = max(e(:));

