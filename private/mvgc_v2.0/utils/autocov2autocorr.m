%% autocov2autocorr
%
% Convert autocovariance to autocorrelation
%
% <matlab:open('autocov2autocorr.m') code>
%
%% Syntax
%
%     R = autocov2autocorr(G)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          a covariance matrix or autocovariance sequence
%
% _output_
%
%     R          a correlation matrix or autocorrelation sequence
%
%% Description
%
% Calculates autocorrelation |R| from autocovariance |G|.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function R = autocov2autocorr(G)

assert(ismatrix(G) || ndims(G) == 3,'autocovariance matrix must have 2 or 3 dimensions');
[n,n1,q1] = size(G);
assert(n1 == n,'autocovariance matrix has bad shape');

R = zeros(n,n,q1);
d = 1./sqrt(diag(G(:,:,1)));
for k = 1:q1
    R(:,:,k) = bsxfun(@times,bsxfun(@times,d,G(:,:,k)),d');
end

%{
D = diag(1./sqrt(diag(G(:,:,1))));
for k = 1:q1
    R(:,:,k) = D*G(:,:,k)*D;
end
%}
