function [ bins, edges ] = discretize_oct(x, n)
%%DISCRETIZE Group numeric data into equally populated bins.
%
%   B = DISCRETIZE(X, N) splits the (rounded) range of X into N equally
%   populated bins and returns the indices of the bins that the elements of X
%   fall into.
%
%   [B, E] = DISCRETIZE(X, N) also returns the edges of the bins.
%
% This function attempts to partially replicate the functionality of Matlab's
% 'discretize' function, but using quantiles instead of even-sized bins.
%
% Pedro Mediano, Jan 2021

e = quantile(x, linspace(0, 1, n + 1));
bins = arrayfun(@(z) find(z >= e, 1, 'last'), x);
bins(bins > n) = n;

if nargout > 1
  edges = e;
end

