%% MVGC Multivariate Granger Causality Toolbox - Release Notes
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% v1.0 2012-10-09
%
% * Initial release
%
%% v2.0 ????-??-??
%
% * gc/autocov_to_smvgc.m and gc/autocov_to_spwcgc.m have been rewritten to use
%   a new, more efficient algorithm for conditional spectral GC calculation, based
%   on a transformed reduced transfer function. This obviates the need for the
%   autocov/CPSD transforms, so the functions core/autocov_xform and
%   core/cpsd_xform have been deprecated.
%
% * gc/autocov_to_pwcgc.m, gc/GCCA_compat/GCCA_tsdata_to_pwcgc.m and
%   gc/subsample/permtest_tsdata_to_pwcgc.m: loop vectorised for efficiency
%   (thanks to A. Erramuzpe Aliaga for the suggestion).
%
% * core/tsdata_to_infocrit.m: new (much faster) OLS algorithm.
%
% * core/tsdata_to_var.m: fixed bug in LWR algorithm that was returning
%   incorrect residuals covariance matrix (this only really affected the
%   single-lag case).
%
% * core/var_to_cpsd.m: enforce Hermitian CPSD (also ensures auto-spectra are
%   real).
%
% * core/tsdata_to_cpsd.m: frequency resolution now defaults to number of
%   observations, window size defaults to Matlab 'pwelsh' default (i.e. 1/8
%   times number of observations).
%
% * core/var_to_autocov.m: check for single precision input if Lyapunov solver
%   fails. Don't balance in eigenvalue calculation. Use 1-norm for relative
%   error check. Now default (automatic calculation of required lags) is
%   specified by empty 'acmaxlags' parameter, and the special case acmaxlags ==
%   'c' just returns the (zero-lag) covariance matrix of the associated VAR(1).
%
% * core/tsdata_to_autocov.m: correction to algorithm for multi-trial case: take
%   mean of per-trial autocovs.
%
% * utils/parcov.m: NEW utility function to compute partial covariance.
%
% * utils/var2itrfun.m: NEW utility function to compute inverse transfer
%   function.
%
% * utils/plot_pw.m, utils/plot_spw.m: warn about complex values, and set them
%   to NaNs.
%
% * utils/bfft.m, utils/bifft.m deprecated (redundant, since Matlab fft, ifft
%   have a 'dim' parameter).
%
% * utils/warn_supp.m, utils/warn_test.m and utils/warn_if.m deprecated because
%   they're confusing and not so useful. All affected code amended accordingly.
%
% * utils/var_specrad.m: don't balance in eigenvalue calculation.
%
% * utils/warn_test.m: code clean-up.
%
% * utils/legacy/randi/randi.m amended to take 3 arguments.
%
% * stats/demean.m: use 'bsxfun' for efficient vectorisation.
%
% * stats/empirical_confints.m: check for too little variation.
%
% * stats/mvgc_confints.m and stats/empirical_confints.m: replaced alpha with
%   alpha/2, since (as the function names imply) we should be computing
%   confidence /intervals/, not bounds!
%
% * startup.m: Warning added not to add full MVGC directory hierarchy to Matlab
%   search path, due to number of users who seem to want to do this. Also added
%   warning about single precision data.
%
% * startup.m: Workaround for Matlab Linux "static TLS" bug
%   (see http://www.mathworks.de/support/bugreports/961964).
%
% * Reduced verbosity of subsampling routines.
%
% * Some general code clean-up.
%
% * Minor amendments to documentation.
%
%%
