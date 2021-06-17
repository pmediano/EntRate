%*********************************************************
%
%
%*********************************************************
if ~exist('OCTAVE_VERSION', 'builtin')
  addpath ..;
end

%% Test with unit-length sequences
assert(LZ76([true]) == 1);
assert(LZ76([false]) == 1);


%% Examples from the Kaspar & Schuster (1987) paper
assert(LZ76([0 0 0 0 0 0 0] < 0.5) == 2);
assert(LZ76([0 1 0 1 0 1 0 1] < 0.5) == 3);


%% Example from the Lempel & Ziv (1976) paper
S = [0 0 0 1 1 0 1 0 0 1 0 0 0 1 0] < 0.5;
assert(LZ76(S) == 6);

