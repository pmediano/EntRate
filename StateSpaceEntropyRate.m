function [ H, band_H ] = StateSpaceEntropyRate(X, Fs, downsampling, band, varmomax)
%%STATESPACEENTROPYRATE LZ-style complexity using a state-space model
%
%     H = STATESPACEENTROPYRATE(X, FS) computes the average entropy rate in X,
%     where X is a matrix of size [D,T,M] containing time series of length T,
%     for D channels and M trials each, and FS is the sampling rate. One model
%     is trained for each channel (using all trials), and entropy rate is
%     averaged across channels. If NDIMS(X) == 2, it is assumed that M = 1.
%
%     H = STATESPACEENTROPYRATE(X, FS, DOWNSAMPLING) specifies whether to
%     downsample the data to a maximum of 200Hz, recommended (can be 'yes' or
%     'no', default: yes).
%
%     [H, BAND_H] = STATESPACEENTROPYRATE(X, FS, DOWNSAMPLING, BAND), where
%     BAND is a K-by-2 array with K frequency bands, returns a length-K array
%     with spectrally decomposed entropy rate over the requested bands (an
%     upper limit of 'inf' defaults to the Nyquist freq).
%
%     [...] = STATESPACEENTROPYRATE(...,VARMOMAX) tries VAR models of order up
%     to VARMOMAX to estimate the optimal SS model order. The default is 20,
%     which works well for typical MEEG data, but not for fMRI. For fMRI, a
%     reasonable heuristic is 8/TR (i.e. VARMOMAX = 2 for standard 2s fMRI).
%
% Example:
%   To compute entropy rate decomposed in the usual frequency bands (delta,
%   theta, etc), use:
%
%   StateSpaceEntropyRate(X, Fs, 'yes', [1, 4; 4, 8; 8, 12; 12, 25]);
%
%   To compute broadband entropy rate in BOLD signals, use:
%
%   StateSpaceEntropyRate(X, 1/TR, [], [], round(8/TR))
%
% This code is based on Lionel Barnett's MVGC toolbox.
%
% Pedro Mediano, Oct 2020

%% Parameter checks and initialisation
if isempty(X) || ~isnumeric(X)
  error('Input must be a 1, 2, or 3D data matrix');
end
if nargin < 2 || isempty(Fs)
  error('Not enough input arguments. You must provide a sampling frequency.');
end
[D, T, M] = size(X);
H = 0;
pmax = 20;

if T < D
  warning(['Your data has more dimensions than timesteps. ', ...
           'Maybe you need to transpose it?']);
end

% Load Octave packages if relevant
if exist('OCTAVE_VERSION', 'builtin') && ~exist('fcdf', 'file')
  pkg('load', 'statistics');
end

% Silently add MVGC to current path
p = mfilename('fullpath');
addpath(strrep(p, 'StateSpaceEntropyRate', 'private/mvgc_v2.0'));
evalc('mvgc_startup;');

if nargin < 3 || isempty(downsampling), downsampling = 'yes'; end
if nargin < 4 || isempty(band),         band         = [];    end
if nargin < 5 || isempty(varmomax),     varmomax     = 20;    end

if nargout < 2 && ~isempty(band)
  error(['Not enough output arguments. When using frequency decomposition, ', ...
        'please use\\[H, band_H] = StateSpaceEntropyRate(...)']);
end

if strcmp(downsampling, 'yes') && Fs > 200
  k = floor(Fs/200);
  Fs = Fs/k;
  X = X(:,1:k:end,:);
end

if ~isempty(band)
  band_H = zeros([size(band,1), 1]);
end


%% Loop over channels and compute entropy rate
for d=1:D
  y = X(d,:,:);

  % Entropy function for a multivariate normal
  % MVGC logdet() is faster, safer and more accurate than log(det())
  H_fun = @(C) 0.5*logdet(2*pi*exp(1)*C);

  % Select VAR model order
  [varmoaic,varmobic,varmohqc,varmolrt] = tsdata_to_varmo(y,varmomax,'LWR',[],false);

  % Select SS model order
  if varmohqc < 1
    % Force the algorithm to fit a SS model even if varmo gives up
    varmohqc = 1;
    ssmo = 2;
  else
    [ssmo,~] = tsdata_to_sssvc(y, 2*varmohqc, [], []);
  end

  % Fit SS model
  % 2*AIC is Bauer's recommendation for past/future window... personally I prefer
  % 2*HQC (HQC = Hannan-Quinn IC, generally sits between AIC and BIC)
  failed = false;
  [A,C,K,V] = tsdata_to_ss(y, 2*varmohqc, ssmo);
  assert(all(size(V) == [1, 1]));  % Check we fit model to 1D data
  try
    info = ss_info(A,C,K,V,0);
    % Manually unset SS error related to negative multi-info (which is often
    % triggered by numerical errors even in 1D models)
    info.error = bitset(info.error, 7, 0);
  catch
    H = nan;
    failed = true;
  end
  if info.error
    H = nan;
    failed = true;
  end

  % Calculate the entropy rate from the SS residuals covariance matrix
  H = H + H_fun(V);


  %% Decompose entropy rate into bands, if requested
  if ~isempty(band)

    if failed
      band_H = nan*ones(size(band_H));
    else

      % Compute CPSD
      fres = 1000;  % Resolution in frequency space
      S = ss_to_cpsd(A,C,K,V,fres);
      H_freq = shiftdim(arrayfun(@(i) H_fun(S(:,:,i)),1:(fres+1)), -1);

      for j=1:size(band,1)
        if isinf(band(j,2))
          band(j,2) = floor(Fs/2.0);
        end

        % Integrate (average) over the requested band (3 is dimension to integrate over)
        band_H(j) = band_H(j) + bandlimit(H_freq, 3, Fs, band(j,:));
      end

    end
  end

end


%% Divide by number of channels and return
H = H/D;
if nargout > 1
  band_H = band_H/D;
end

