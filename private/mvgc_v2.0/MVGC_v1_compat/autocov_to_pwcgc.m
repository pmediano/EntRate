%% autocov_to_pwcgc
%
% Calculate pairwise-conditional time-domain MVGCs (multivariate Granger causalities)
%
% <matlab:open('autocov_to_pwcgc.m') code>
%
%% Syntax
%
%     F = autocov_to_pwcgc(G)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%
% _output_
%
%     F          Granger causality matrix
%
%% Description
%
% Returns the  matrix |F| of pairwise-conditional time-domain MVGCs
%
% <<eq_mvgc_pwc.png>>
%
% (where |[ij]| denotes omission of the |ij|-th variables) between all
% pairs of variables [[ii_inej.png]] represented in |G|, for a stationary VAR
% process with autocovariance sequence |G|. Note that the first index |i| of
% |F| is the target (causee) variable, the second |j| the source (causal)
% variable. See ref. [1] for details.
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(F,false)|.
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
% <autocov_to_mvgc.html |autocov_to_mvgc|> |
% <autocov_to_var.html |autocov_to_var|> |
% <mvgc_demo.html |mvgc_demo|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function F = autocov_to_pwcgc(G)

n = size(G,1);
F = nan(n);

% full regression

wstate = warning('off','all'); lastwarn('');
[~,SIG] = autocov_to_var(G);
wmsg = lastwarn; warning(wstate);
if ~isempty(wmsg), fprintf(2,'WARNING in full regression (check output of ''var_info''): %s\n',wmsg); end
if isbad(SIG),     fprintf(2,'ERROR in full regression (check output of ''var_info'')\n'); return; end % show-stopper!

LSIG = log(diag(SIG));

for j = 1:n

    % reduced regression
    
    jo = [1:j-1 j+1:n]; % omit j

    wstate = warning('off','all'); lastwarn('');
    [~,SIGj] = autocov_to_var(G(jo,jo,:));
    wmsg = lastwarn; warning(wstate);
    if ~isempty(wmsg), fprintf(2,'WARNING in reduced regression for source node %d (check output of ''var_info''): %s\n',j,wmsg); end
    if isbad(SIGj),    fprintf(2,'ERROR in reduced regression for source node %d (check output of ''var_info'')\n',j); continue; end % show-stopper!
    
    F(jo,j) = log(diag(SIGj))-LSIG(jo); % v2.0 - loop vectorised (thanks to A. Erramuzpe Aliaga)
    %{
    LSIGj = log(diag(SIGj));
    for ii=1:n-1;
        i = jo(ii);
        F(i,j) = LSIGj(ii)-LSIG(i);
    end
    %}
end
