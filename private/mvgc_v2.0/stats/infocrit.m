%% infocrit
%
% Calculate Akaike, Bayesian and Hannan-Quinn information criteria
%
% <matlab:open('infocrit.m') code>
%
%% Syntax
%
%     [aic,bic] = infocrit(L,k,m)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     L          log-likelihood
%     k          number of free parameters
%     m          number of observations
%     HandT      use Hurvich and Tsai AIC small-sample bias correction
%
% _output_
%
%     aic        AIC value
%     bic        BIC value
%     hqc        HQC value
%
%% Description
%
% Calculates the Akaike (AIC), Bayesian (BIC) and Hannan-Quinn (HQC) criteria
% from the log-likelihood of a model [1]. |L|, |k| and |m| may be scalars, or
% matching vectors or matrices. For the Akaike information criterion, the
% small-sample bias-corrected form ("AICc") of Hurvich and Tsai [2] is
% calculated.
%
%% References
%
% [1] K. P. Burnham and D. R. Anderson, "Model Selection and Multimodel
% Inference: A Practical Information-Theoretic Approach", _2nd ed._,
% Springer-Verlag, 2002.
%
% [2] C. M. Hurvich and C.-L. Tsai, "Regression and time series model selection
% in small samples", _Biometrika_, 76, 1989.
%
%% See also
%
% <tsdata_to_infocrit.html |tsdata_to_infocrit|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [aic,bic,hqc] = infocrit(L,k,m,HandT)

if nargin <  4 || isempty(HandT), HandT = false; end

K = k./m;
if HandT
    fac = m./(m-k-1);
    aic = -2*L + 2*K.*fac; % Akaike's AIC with Hurvich and Tsai's small-sample correction
    aic(fac <= 0) = NaN;
else
    aic = -2*L + 2*K;               % Akaike's AIC
end
bic = -2*L + K.*log(m);             % Schwarz' BIC
hqc = -2*L + 2*K.*log(log(m));      % Hannan-Quinn IC
