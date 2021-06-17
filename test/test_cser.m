%*********************************************************
%
%
%*********************************************************
if ~exist('OCTAVE_VERSION', 'builtin')
  addpath ..;
end
N = 10000;
nb_repeats = 10;
tol = 1e-2;

H_fun = @(S) 0.5*log(det(2*pi*exp(1)*S));

%% IID data
S_vec = [0.1, 0.5, 1, 2, 5];
for S=S_vec
  H_tmp = zeros([1 nb_repeats]);
  for i=1:nb_repeats
    H_tmp(i) = StateSpaceEntropyRate(sqrt(S)*randn([1, N]), 1);
  end
  H = mean(H_tmp);
  assert(abs(H - H_fun(S)) < tol);
end


%% Multi-trial data
S_vec = [0.1, 0.5, 1, 2, 5];
for S=S_vec
  H = StateSpaceEntropyRate(sqrt(S)*randn([1, N, nb_repeats]), 1);
  assert(abs(H - H_fun(S)) < tol);
end


%% Multi-channel data
H_tmp = zeros([1 nb_repeats]);
for i=1:nb_repeats
  H_tmp(i) = StateSpaceEntropyRate([randn([1, N]); 2*randn([1, N])], 1);
end
H = mean(H_tmp);
assert(abs(H - (H_fun(1)+H_fun(4))/2) < tol);


%% Spectral decomposition of white noise
% In white noise, all equally-sized bands should have the same CSER
[H, bH] = StateSpaceEntropyRate(randn([1, N, nb_repeats]), 1, [], ...
                                [0.05, 0.10; 0.12, 0.17; 0.20, 0.25]);
assert(all(abs(bH(1) - bH) < tol));
