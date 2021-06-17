%% autocov_to_spwcgc
%
% Calculate pairwise-conditional frequency-domain MVGCs (spectral multivariate Granger causalites)
%
% <matlab:open('autocov_to_spwcgc.m') code>
%
%% Syntax
%
%     [f,fres] = autocov_to_spwcgc(G,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     fres       frequency resolution
%
% _output_
%
%     f          spectral Granger causality matrix
%
%% Description
%
% Returns the  matrix |f| of pairwise-conditional frequency-domain
% (spectral) MVGCs
%
% <<eq_smvgc_pwc.png>>
%
% (where |[ij]| denotes omission of the |ij|-th variables) between all
% pairs of variables [[ii_inej.png]] represented in |G|, for a stationary VAR
% process with autocovariance sequence |G|. The first index |i| of
% |f| is the target (causee) variable, the second |j| the source (causal)
% variable and the third indexes the frequency. See ref. [1] for details.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. If |fres| is not supplied it is calculated optimally
% as the number of autocovariance lags. Call |freqs =
% <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the sampling
% rate, to get a corresponding vector |freqs| of frequencies on |[0,fs/2]|.
%
% A new algorithm is used based on transformed reduced transfer functions.
% Separate estimation steps for the reduced regressions (which is known to be
% problematic [2,*]) are not required, resulting in good efficiency and
% accuracy.
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(f,false)|.
%
% For details of the algorithm, see <autocov_to_smvgc.html |autocov_to_smvgc|> and [1].
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
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_pwcgc.html |autocov_to_pwcgc|> |
% <autocov_to_var.html |autocov_to_var|> |
% <var2trfun.html |var2trfun|> |
% <parcov.html |parcov|> |
% <sfreqs.html |sfreqs|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [f,fres] = autocov_to_spwcgc(G,fres,~)

% v2.0 - completely rewritten using a new algorithm; doesn't need
% 'autocov_xform()' anymore (so old 'useFFT' parameter now redundant, replaced
% with '~' for backward compatibility). See autocov_to_smvgc.m.

[n,~,q1] = size(G);
if nargin < 2 || isempty(fres);
    fres = q1;
end

h = fres+1;
f = nan(n,n,h);

wstate = warning('off','all'); lastwarn('');
[A,SIG] = autocov_to_var(G); % full regression
wmsg = lastwarn; warning(wstate);
if ~isempty(wmsg), fprintf(2,'WARNING in full regression (check output of ''var_info''): %s\n',wmsg); end
if isbad(A),       fprintf(2,'ERROR in full regression (check output of ''var_info'')\n'); return; end % show-stopper!

wstate = warning('off','all'); lastwarn('');
H = var2trfun(A,fres); % full transfer function
wmsg = lastwarn; warning(wstate);
if ~isempty(wmsg), fprintf(2,'WARNING in transfer function calculation: %s\n',wmsg); end
if isbad(H),       fprintf(2,'ERROR in transfer function calculation\n'); return; end % show-stopper!

for j = 1:n
    oj = [1:j-1 j+1:n]; % omit j

    wstate = warning('off','all'); lastwarn('');
    [AR,SIGR] = autocov_to_var(G(oj,oj,:));  % reduced regression
    wmsg = lastwarn; warning(wstate);
    if ~isempty(wmsg), fprintf(2,'WARNING in reduced regression for source node %d (check output of ''var_info''): %s\n',j,wmsg); end
    if isbad(AR),      fprintf(2,'ERROR in reduced regression for source node %d (check output of ''var_info'')\n',j); continue; end % show-stopper!

    BR = var2itrfun(AR,fres); % reduced inverse transfer function

    for ii = 1:n-1;
        i  = oj(ii);         % i index in omitted j indices
        oi = [1:i-1 i+1:n];  % omit i

        SRii  = SIGR(ii,ii); % reduced spectrum is flat!
        LSRii = log(SRii);

        PSIGSR = chol(parcov(SIG,oi,i),'lower'); % partial covariance square root

        for k = 1:h
            HRk = BR(ii,:,k)*H(oj,oi,k)*PSIGSR; % transformed reduced transfer function
            f(i,j,k) = LSRii - log(SRii-HRk*HRk');
        end
    end
end
