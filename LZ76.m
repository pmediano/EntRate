%%LZ76 Lempel-Ziv complexity, 1976 version.
%
%    C = LZ76(X) computes the unnormalised LZ complexity of the 1D (row or
%    column) boolean array X. The recommended way to normalise LZ complexity
%    is C*log2(length(X))/length(X), to obtain an estimate of the entropy rate
%    of X. Surrogate normalisation and Hilbert transform are discouraged.
%
% NOTE: this function needs to be mex'd before use. With a functional mex
% installation, run (from the Matlab prompt):
%
%    mex COPTIMFLAGS="-O3" LZ76.c
%
% Reference:
%   Kaspar and Schuster (1987). Easily-calculable measure for the complexity of
%   spatiotemporal patterns.
%   
% Pedro Mediano and Andrea Luppi, Oct 2020
