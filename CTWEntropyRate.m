function [ H ] = CTWEntropyRate(X, symbols, vmo)
%%CTWENTROPYRATE LZ-style complexity using a context tree-weighted algorithm
%
%     H = CTWENTROPYRATE(X) discretises the signal X and computes its entropy
%     rate using a CTW Variable-order Markov Model (VMM), where X is a matrix
%     of size [D,T,M] containing time series of length T, for D channels and M
%     trials each. One model is trained for each channel (using all trials),
%     and entropy rate is averaged across channels. If NDIMS(X) == 2, it is
%     assumed that M = 1.
%
%     H = CTWENTROPYRATE(X, S) discretises the signal into S different discrete
%     symbols before feeding it into CTW. A large number of symbols may offer
%     greater resolution, but can also worsen statistics and require more data
%     to converge (default: 2).
%
%     H = CTWENTROPYRATE(X, S, VMO) uses a maximum VMM model order of VMO
%     (default: 100).
%
% This code is based on Ron Begleiter's VMM CTW library.
%
% Pedro Mediano, Dec 2020

%% Parameter checks and initialisation
if isempty(X) || ~(isnumeric(X) || islogical(X))
  error('Input must be a 1, 2, or 3D data matrix');
end
[D, T, M] = size(X);
H = 0;
if T < D
  warning(['Your data has more dimensions than timesteps. ', ...
           'Maybe you need to transpose it?']);
end

% Silently add the VMM jars and base functions to current path
p = strrep(mfilename('fullpath'), 'CTWEntropyRate', 'private/');
javaaddpath([p, 'vmm/vmm.jar'])
javaaddpath([p, 'vmm/trove.jar'])
if exist([p, '../../elph_base'], 'dir')
  addpath([p, '../../elph_base']);
end

if nargin < 2 || isempty(symbols), symbols =   2; end
if nargin < 3 || isempty(vmo),     vmo     = 100;    end

alphabet_size = symbols + 1;


%% Find whether data needs to be discretized
% If data is already discrete and has the desired number of symbols (or less),
% then do not discretize further
if isdiscrete(X) && length(unique(X(:))) <= symbols
  disc_fun = @(z) z;

else
  % Otherwise, discretize in quantiles
  disc_fun = @(z) discretize_oct(z, symbols);
end


%% Loop over channels and compute entropy rate
for d=1:D

  % Initialise VMM-CTW model
  mdl = javaObject('vmm.algs.DCTWPredictor');
  mdl.init(alphabet_size, vmo);

  % Make cell array with trials, discretise them and make them java strings
  trials  = num2cell(X(d,:,:), 2);
  dtrials = cellfun(@(y) javaObject('java.lang.String', chararray(disc_fun(y))), ...
                    trials, 'UniformOutput', false);

  % Feed each trial separately into the CTW training algorithm
  cellfun(@(t) mdl.learn(t), dtrials, 'UniformOutput', false);

  % Evaluate their log-likelihood and average to get entropy rate
  H = H + mean(cellfun(@(t) mdl.logEval(t), dtrials))/(T*D);

end

end

function c = chararray(s)
% Auxiliary function to tokenise a string or int array and remap it to
% character codes starting from 1. To print string, use sprintf('%d', c).
% Adapted from Ron Begleiter's @alphabet class.

ab_str = sort(unique(s));
c = char();
for i=1:length(ab_str)
  c(find(s==ab_str(i))) = i; % NOTE: index of first char = "1" (in JAVA it should be "0") 
                             % thus, size(ab) = num of symbols+1 
end

end

