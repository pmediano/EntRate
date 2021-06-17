%% var_to_autocov
%
% Return autocovariance sequence for a VAR model
%
% <matlab:open('var_to_autocov.m') code>
%
%% Syntax
%
%     [G,info] = var_to_autocov(A,V,acmaxlags,acdectol,aitr,maxiters,maxrelerr)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     V          residuals covariance matrix
%     acmaxlags  maximum autocovariance lags to calculate, empty for automatic calculation (default), or
%                (special case) 'c' for the covariance matrix of the process
%     acdectol   autocovariance decay tolerance (default: 1e-8)
%     aitr       use "accelerated iterative" Lyapunov solver algorithm, else (default) a Schur algorithm (see Description)
%     maxiters   maximum iterations if using iterative algorithm; see 'dlyap_aitr' for defaults
%     maxrelerr  maximum relative error if using iterative algorithm; see 'dlyap_aitr' for defaults
%
% _output_
%
%     G          autocovariance sequence
%     info       info structure, with fields (some may not be present):
%         error      error number (0 if no error)
%         errmsg     error message string
%         warnings   number of warnings (0 if no warnings)
%         warnmsg    warning mesage strings (cell array)
%         rho        VAR spectral radius
%         iters      number of iterations if using iterative algorithm, else 0
%         acrelerr   relative error of associated VAR(1) solution
%         acminlags  minimum lags required to achieve specified autocovariance decay factor
%         aclags     actual number of autocovariance lags calculated
%
%% Description
%
% Returns autocovariance sequence |G| defined as [[ii_acseq.png]]
% for a VAR model with coefficients |A| and (positive-definite) residual
% covariance matrix |V|, by "reverse-solving" the Yule-Walker equations
%
% <<eq_yweqs.png>>
%
% (where  [[ii_Sigma.png]] = |V|). The algorithm solves for the associated
% VAR(1) covariance matrix - a discrete-time Lyapunov equation [2,3] - and then
% calculates higher lags recursively [1].
%
% Errors, warnings and diagnostics are returned in the |info| struct. Use the
% routine <var_info.html |var_info|> (which should _always_ be called after this
% function) to display this information. Possible errors are
%
%     info.error     info.errmsg
%     ---------------------------------------------------------------------------
%         0          (no error, no message)
%         1          unstable VAR
%         2          residuals covariance matrix not positive-definite
%         3          Lyapunov equation solver failed for some reason
%         4          covariance matrix of associated VAR(1) not positive-definite
%     ---------------------------------------------------------------------------
%
% For a stable VAR the the spectral radius |info.rho|  (see <specnorm.html
% |specnorm|>) must be < 1; this may be considered a crude unit root test for
% stationarity [1]. Then the autocovariance sequence decays approximately
% exponentially, by a factor equal to |info.rho|. The minimum number of lags
% required to achieve the specified autocovariance decay tolerance |acdectol| is
% calculated as |info.acminlags|, so that |info.rho^info.acminlags < acdectol|.
% The actual number of lags |info.aclags| to which autocovariance is calculated
% is then set to the minimum of |info.acminlags| and the specified maximum
% number of lags, |acmaxlags| (if |acmaxlags| is not supplied - the recommended
% option - it defaults to |info.acminlags|). A warning is issued if |info.aclags
% < info.acminlags|. In this case there is no guarantee that MVGCs -
% particularly in the spectral domain - will be accurate. However, if the
% spectral radius of the VAR model is close to 1, so that |info.acminlags| is
% unfeasibly large, there may be no alternative [note that most Granger
% causality libraries effectively set |info.aclags| to the model order]. The
% special case |acmaxlags == 'c'| just returns the (zero-lag) covariance matrix
% for the associated VAR(1) model [1].
%
% The covariance matrix for the associated VAR(1) model is checked for
% positive-definitiveness. If this check fails, it may be an indication of an
% ill-conditioned VAR (possibly because residuals variances are too small, there
% is excessive co-linearity in the data and/or the process is borderline
% stationary); we stress that the |info.error| field _must_ be checked by the
% caller: 0 signifies success, > 0 signifies an error with corresponding message
% in |info.errmsg|. Attention should also be paid to any warnings in
% |info.WARNINGn|.
%
% The |aitr| flag specifies a fast (but experimental) "accelerated" iterative
% Lyapunov solver algorithm (see <dlyap_aitr.html |dlyap_aitr|>). If not set
% (default), the Control System Toolbox <matlab:doc('dlyap') |dlyap|> solver
% routine is used if available, else a roughly equivalent (but slower) scripted
% algorithm based on Schur decomposition (see <matlab:open('dlyap.m') |dlyap.m|>
% in the |utils/control| directory).
%
% *_Note:_* This routine may suffer from numerical instability if the input is
% single precision floating point. Please ensure that your input is double
% precision.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] X. Kitagawa, An algorithm for solving the matrix equation X = F X F'+S,
%     _Internat. J. Control_ 25(5), 1977.
%
% [3] S. Hammarling, "Numerical solution of the stable, non-negative definite
%     Lyapunov equation", _IMA J. Numer. Anal._ 2, 1982.
%
%% See also
%
% <specnorm.html |specnorm|> |
% <var_info.html |var_info|> |
% <matlab:doc('dlyap') |dlyap (ControlSystem Toolbox)|> |
% <dlyap.html |dlyap (scripted Schur algorithm)|> |
% <dlyap_aitr.html |dlyap_aitr|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [G,info] = var_to_autocov(A,V,acmaxlags,acdectol,aitr,maxiters,maxrelerr)

%global have_dlyap;

% default parameters

if nargin < 3,                       acmaxlags = [];    end % v2.0 - calculate maximum lags automatically
if nargin < 4 || isempty(acdectol),  acdectol  = 1e-8;  end % autocovariance decay tolerance
if nargin < 5 || isempty(aitr),      aitr      = false; end % use "accelerated" iterative Lyapunov equation solver

% iterative algorithm only: ensure defaults for utils/dlyap_aitr.m.

if nargin < 6, maxiters  = []; end
if nargin < 7, maxrelerr = []; end

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn = p*n;
pn1 = (p-1)*n;

[nn1,nn2] = size(V);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

% initialise info struct
info.error     = 0;
info.errmsg    = '';
info.warnings  = 0;
info.warnmsg   = cell(0,1);
info.rho       = NaN;
info.mii       = NaN;
info.mmii      = NaN;
info.iters     = NaN;
info.acrelerr  = NaN;
info.acminlags = NaN;
info.aclags    = NaN;

G = [];

% construct VAR coefficients for associated VAR(1)

A1 = [reshape(A,n,p*n); eye(pn1) zeros(pn1,n)]; % "companion matrix" for A

% calculate spectral radius

info.rho = max(abs(eig(A1,'nobalance'))); % v2.0 - don't balance!

if info.rho > 1-eps % v2.0 - to be on the safe side
    info.error = 1;
    info.errmsg = 'unstable VAR';
    return
end

% construct residual covariances for associated VAR(1)

info.mii  = multiinfo(V);     % multi-information
info.mmii = multiinfo(n,true);  % multi-information for uniform random n x n correlation matrix, for comparison

if maxabs(triu(V,1)-triu(V',1)) > eps || ~isposdef(V);
    info.error = 2;
    info.errmsg = 'residuals covariance matrix not symmetric, positive definite';
    return;
end

V1 = [V zeros(n,pn1); zeros(pn1,n) zeros(pn1)];

% solve the Lyapunov equation for the covariance matrix of the associated VAR(1)

try
    if aitr
        [G1,info.iters] = dlyap_aitr(A1,V1,maxiters,maxrelerr); % experimental: fast, but needs more testing
    else
%       G1 = dlyap(A1,V1);           % dlyap seems to work better here without balancing, which seems to break positive-definitiveness
        G1 = lyapslv('D',A1,[],-V1); % sometimes. However lyapslv is not an official interface, so this could conceivably break in future.
    end
catch except
    info.error = 3;
    if isa(A,'single') || isa(V,'single') % v2.0 single precision data may cause this error
        info.errmsg = sprintf('Lyapunov equation solver failed: %s\n    This may be due to single precision input. Please convert your data to double precision!', except.message);
    else
        info.errmsg = ['Lyapunov equation solver failed: ' except.message];
    end
    return
end

info.acrelerr = norm(A1*G1*A1'-G1+V1,1)/(1+norm(G1,1)+norm(V1,1)); % this should be small (see below) v2.0 use 1-norm

maxacrelerr = 1e-8; % probably shouldn't be hard-coded :-/
if info.acrelerr > maxacrelerr
    info.warnings = info.warnings+1;
    info.warnmsg{info.warnings} = sprintf('large relative error = %g (tolerance = %g)',info.acrelerr,maxacrelerr);
end

% estimate number of autocov lags

info.acminlags = ceil(log(acdectol)/log(info.rho)); % minimum lags to achieve specified tolerance

% set actual number of autocov lags

if isempty(acmaxlags) % use minimum acceptable lags (recommended)
   info.aclags = info.acminlags;
else
    assert(isscalar(acmaxlags),'bad ''acmaxlags'' parameter')
    if ischar(acmaxlags) && lower(acmaxlags) == 'c' % special case - covariance matrix of associated VAR(1)
        G = G1;
        return;
    end
    assert(isint(acmaxlags),'bad ''acmaxlags'' parameter (must be an integer)')
    if     acmaxlags <= 0  % use exactly -acmaxlags lags, and don't issue a warning (not encouraged, hence undocumented!)
        info.aclags = -acmaxlags;
    else                   % use at most acmaxlags lags
        info.aclags  = min(info.acminlags,acmaxlags);
        if info.aclags < info.acminlags
            info.warnings = info.warnings+1;
            info.warnmsg{info.warnings} = 'too few autocovariance lags';
        end
    end
end

if ~isposdef(G1);
    info.error = 4;
    info.errmsg = 'associated VAR(1) covariance matrix not positive-definite';
    return
end

q = info.aclags;
q1 = q+1;

% we already have p-1 lags; if that's enough, truncate (if necessary) and return

if q < p
    G = reshape(G1(1:n,1:n*q1),n,n,q1);
    return
end

% calculate recursively from associated VAR(1) solution, from p lags up to q lags

G = cat(3,reshape(G1(1:n,:),n,n,p),zeros(n,n,q1-p));    % autocov forward  sequence
B = [zeros((q1-p)*n,n); G1(:,end-n+1:end)];             % autocov backward sequence
A = reshape(A,n,pn);                                    % coefficients
for k = p:q
    r = q1-k;
    G(:,:,k+1) = A*B(r*n+1:r*n+pn,:);
    B((r-1)*n+1:r*n,:) = G(:,:,k+1);
end
