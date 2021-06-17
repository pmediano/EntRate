%% autocov_to_mvgc
%
% Calculate conditional time-domain MVGC (multivariate Granger causality)
%
% <matlab:open('autocov_to_mvgc.m') code>
%
%% Syntax
%
%     F = autocov_to_mvgc(G,x,y)
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
%
% _output_
%
%     F          Granger causality
%
%% Description
%
% Returns the time-domain MVGC
%
% <<eq_mvgc.png>>
%
% from the variable |Y| (specified by the vector of indices |y|) to the
% variable |X| (specified by the vector of indices |x|), conditional on all
% other variables |Z| represented in |G|, for a stationary VAR process with
% autocovariance sequence |G|. See ref. [1] for details.
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
% <autocov_to_var.html |autocov_to_var|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function F = autocov_to_mvgc(G,x,y)

n = size(G,1);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)
r = [x z];

F = NaN;

% full regression

wstate = warning('off','all'); lastwarn('');
[~,V] = autocov_to_var(G);
wmsg = lastwarn; warning(wstate);
if ~isempty(wmsg), fprintf(2,'WARNING in full regression (check output of ''var_info''): %s\n',wmsg); end
if isbad(V),     fprintf(2,'ERROR in full regression (check output of ''var_info'')\n'); return; end % show-stopper!

% reduced regression

wstate = warning('off','all'); lastwarn('');
[~,VR] = autocov_to_var(G(r,r,:));    % reduced regression
wmsg = lastwarn; warning(wstate);
if ~isempty(wmsg), fprintf(2,'WARNING in reduced regression (check output of ''var_info''): %s\n',wmsg); end
if isbad(VR),    fprintf(2,'ERROR in reduced regression (check output of ''var_info'')\n'); return; end % show-stopper!

xr = 1:length(x);
F = logdet(VR(xr,xr))-logdet(V(x,x));
