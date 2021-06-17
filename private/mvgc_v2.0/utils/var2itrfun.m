%% var2itrfun
%
% Calculate VAR inverse transfer function from VAR coefficients
%
% <matlab:open('var2itrfun.m') code>
%
%% Syntax
%
%     J = var2itrfun(A,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     fres       frequency resolution
%
% _output_
%
%     J          VAR inverse transfer function matrix
%
%% Description
%
% Return inverse transfer function |J| for VAR with coefficients |A|. |fres|
% specifies the frequency resolution. Call |freqs = <sfreqs.html
% sfreqs>(fres,fs)|, where |fs| is the sampling rate, to get a corresponding
% vector |freqs| of frequencies on |[0,fs/2]|.
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
% <var2trfun.html |var2trfun|> |
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_spwcgc.html |autocov_to_spwcgc|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function J = var2itrfun(A,fres)

n = size(A,1);
J = fft(cat(3,eye(n),-A),2*fres,3); % over [0,2*pi)
J = J(:,:,1:(fres+1));              % over [0,pi] only
