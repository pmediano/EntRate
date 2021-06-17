%% autocov_to_smvgc
%
% Calculate conditional frequency-domain MVGC (spectral multivariate Granger causality)
%
% <matlab:open('autocov_to_smvgc.m') code>
%
%% Syntax
%
%     [f,fres] = autocov_to_smvgc(G,x,y,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     fres       frequency resolution (default: automatic)
%
% _output_
%
%     f          spectral Granger causality
%
%% Description
%
% Returns the frequency-domain (spectral) MVGC
%
% <<eq_smvgc.png>>
%
% from the variable |Y| (specified by the vector of indices |y|) to the
% variable |X| (specified by the vector of indices |x|), conditional on all
% other variables |Z| represented in |G|, for a stationary VAR process with
% autocovariance sequence G.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. If |fres| is not supplied it is calculated optimally
% as the number of autocovariance lags. Call |freqs =
% <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the sampling
% rate, to get a corresponding vector |freqs| of frequencies on |[0,fs/2]|.
%
% In the conditional case, a new algorithm is used based on a transformed
% reduced transfer function. A separate estimation step for the reduced
% regression (which is known to be problematic [2,*]) is not required, resulting
% in good efficiency and accuracy.
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(f,false)|.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] Y. Chen, S. L. Bressler and M. Ding, "Frequency decomposition of
% conditional Granger causality and application to multivariate neural
% field potential data", _J. Neurosci. Methods_, 150, 2006.
%
% [*] In our experience the "partition matrix" method in ref. [2] appears to be
% unsound, producing inaccurate results; hence we do not use it here.
%
%% See also
%
% <autocov_to_var.html |autocov_to_var|> |
% <var_to_cpsd.html |var_to_cpsd|> |
% <var2trfun.html |var2trfun|> |
% <parcov.html |parcov|> |
% <sfreqs.html |sfreqs|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [f,fres] = autocov_to_smvgc(G,x,y,fres,~)

% v2.0 - completely rewritten using a new algorithm; doesn't need
% 'autocov_xform()' anymore (so old 'useFFT' parameter now redundant, replaced
% with '~' for backward compatibility)

[n,~,q1] = size(G);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)

if nargin < 4 || isempty(fres), fres = q1; end

h = fres+1;
f = nan(1,h);

wstate = warning('off','all'); lastwarn('');
[A,SIG] = autocov_to_var(G); % full regression
wmsg = lastwarn; warning(wstate);
if ~isempty(wmsg), fprintf(2,'WARNING in full regression (check output of ''var_info''): %s\n',wmsg); end
if isbad(A),       fprintf(2,'ERROR in full regression (check output of ''var_info'')\n'); return; end % show-stopper!

if isempty(z) % unconditional

    wstate = warning('off','all'); lastwarn('');
    [S,H] = var_to_cpsd(A,SIG,fres); % full CPSD and transfer function
    wmsg = lastwarn; warning(wstate);
    if ~isempty(wmsg), fprintf(2,'WARNING in CPSD calculation: %s\n',wmsg); end
    if isbad(S),       fprintf(2,'ERROR in CPSD calculation (covariance matrix not positive definite?\n'); return; end % show-stopper!

    PSIGSR = chol(parcov(SIG,y,x),'lower'); % partial covariance square root

    for k = 1:h
        Hk = H(x,y,k)*PSIGSR; % transformed reduced transfer function
        f(k) = logdet(S(x,x,k)) - logdet(S(x,x,k)-Hk*Hk');
    end

else % conditional

    xz = [x z];
    yz = [y z];
    xr = 1:length(x);

    H = var2trfun(A,fres); % full transfer function

    wstate = warning('off','all'); lastwarn('');
    [AR,SIGR] = autocov_to_var(G(xz,xz,:)); % reduced regression
    wmsg = lastwarn; warning(wstate);
    if ~isempty(wmsg), fprintf(2,'WARNING in reduced regression (check output of ''var_info''): %s\n',wmsg); end
    if isbad(AR),      fprintf(2,'ERROR in reduced regression (check output of ''var_info'')\n'); return; end % show-stopper!

    BR = var2itrfun(AR,fres); % reduced inverse transfer function

    SRxx = SIGR(xr,xr);       % reduced spectrum is flat!
    LDSRxx = logdet(SRxx);

    PSIGSR = chol(parcov(SIG,yz,x),'lower'); % partial covariance square root

    for k = 1:h
        HRk = BR(xr,:,k)*H(xz,yz,k)*PSIGSR; % transformed reduced transfer function
        f(k) = LDSRxx - logdet(SRxx-HRk*HRk');
    end

end
