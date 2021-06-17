%% relerr
%
% Calculate relative error between two arrays
%
% <matlab:open('maxabs.m') code>
%
%% Syntax
%
%     e = relerr(X,Y,)
%
%% Arguments
%
% _input_
%
%     X          an array
%     Y          an array (same dimensions as X)
%     p          norm type, as in Matlab 'norm' function (default: 1-norm)
%
% _output_
%
%     e          the relative error between X and Y
%
%% Description
%
% Returns the relative error |e| between |X| and |Y| (vectorised), using the
% |p|-norm. |NaN| s are *not* ignored.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function e = relerr(X,Y,p)

assert(isequal(size(X),size(Y)),'Arrays must have the same dimensions');

if nargin < 3 || isempty(p), p = 1; end % 1-norm

e = norm(X(:)-Y(:),p)/(1+norm(X(:),p)+norm(Y(:),p));
