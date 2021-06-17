%*********************************************************
%
%
%*********************************************************
if ~exist('OCTAVE_VERSION', 'builtin')
  addpath ..;
end
N = 10000;
nb_repeats = 10;
tol = 3e-2;

% Entropy of a binary Bernoulli variable -- for ground truth
h_bin = @(p) - p*log2(p) - (1-p)*log2(1-p);

%% Fully predictable sequence
% (Adding one different symbol to avoid numerical problems)
H_tmp = zeros([1 nb_repeats]);
for i=1:nb_repeats
  H_tmp(i) = CTWEntropyRate([1 repmat(0, 1, N-1)]);
end
H = mean(H_tmp);
assert(abs(H) < tol);


%% Fully random sequence
H_tmp = zeros([1 nb_repeats]);
for i=1:nb_repeats
  H_tmp(i) = CTWEntropyRate(1*(rand([1, N]) < 0.5));
end
H = mean(H_tmp);
assert(abs(H - 1) < tol, ...
  sprintf('Fully random sequence. Expected 1 bits but got %f.', H));


%% Zero-th order sequence
p = 0.2;
H_tmp = zeros([1 nb_repeats]);
for i=1:nb_repeats
  X = 1*(rand([1, N]) < p);
  H_tmp(i) = CTWEntropyRate(X);
end
H = mean(H_tmp);
assert(abs(H - h_bin(p)) < 5*tol, ...
  sprintf('Zero-th order sequence. Expected %f bits but got %f.', h_bin(p), H));


%% Finite order XOR
p = 0.2;
H_tmp = zeros([1 nb_repeats]);
for i=1:nb_repeats
  X = (rand([1, N]) < p);
  for t=3:N
    X(t) = xor(xor(X(t-2), X(t-1)), X(t));
  end
  H_tmp(i) = CTWEntropyRate(1*X);
end
H = mean(H_tmp);
assert(abs(H - h_bin(p)) < 5*tol, ...
  sprintf('Finite-order XOR. Expected %f bits but got %f.', h_bin(p), H));


%% Finite model order
H = CTWEntropyRate(repmat([0 0 1 1], 1, ceil(N/4)), 2, 2);
assert(abs(H) < tol);

H = CTWEntropyRate(repmat([0 0 1 1], 1, ceil(N/4)), 2, 1);
assert(abs(H - 1) < tol);


%% Invariance to symbol permutation
X = rand([1, N]) < 0.5;
H0 = CTWEntropyRate(1*X);
H1 = CTWEntropyRate(1*(~X));
assert(abs(H0 - H1) < 1e-10);


%% Continuous data: fully random
for k=2:5
  H_tmp = zeros([1 nb_repeats]);
  for i=1:nb_repeats
    H_tmp(i) = CTWEntropyRate(randn([1, N]), k);
  end
  H = mean(H_tmp);
  assert(abs(H - log2(k)) < 2*tol);
end


%% Continuous data: shift invariance
X = zeros([1, N]);
for t=2:N X(t) = 0.8*X(t-1) + randn; end
for k=2:5
  H0 = CTWEntropyRate(X, k);
  H1 = CTWEntropyRate(X + 100, k);
  assert(abs(H0 - H1) < 1e-10);
end

