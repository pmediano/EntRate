%% specnorm
%
% Calculate VAR or VMA spectral norm
%
% <matlab:open('specnorm.m') code>
%
%% Syntax (VAR)
%
%     rho     = specnorm(A)
%     [A,rho] = specnorm(A,newrho)
%
%% Syntax (VMA)
%
%     rho     = specnorm(-B)
%     [B,rho] = -specnorm(-B,newrho)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR (or -VMA) coefficients array - see _NOTE_ below
%     newrho     new VAR spectral norm
%
% _output_
%
%     A          VAR coefficients array
%     rho        VAR spectral norm
%
%% Description
%
% For a VAR we assume a transfer function A(z)^{-1} with
% A(z) = I - A(:,:,1) z - A(:,:,2) z^2 + ..., while for a VMA we assume
% a transfer function B(z) with B(x) = I + B(:,:,1) z + B(:,:,2) z^2 + ...
% For a VAR, supply A, but for a VMA supply -B. For a VAR the process
% is _stable_ iff |rho < 1|, for a VMA the process is _minimum phase_
% iff |rho < 1|.
%
% _First form:_ return the spectral norm |rho| for a VAR/VMA process with
% coefficients matrix |A|. May be used for unit root test (i.e. need |rho < 1|).
%
% _Second form:_, a new value |newrho| for the spectral norm is supplied
% and the coefficients |A| are "decayed" (see function <var_decay.html
% |var_decay|>) so that the new value of the spectral norm becomes
% |newrho|. The new coefficients and old value of the spectral norm are
% returned. Note: |newrho| may be negative, but |newrho < 1| is required
% for the new coefficients to specify a stable process.
%
%% See also
%
% <var_decay.html |var_decay|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [out1,out2] = specnorm(A,newrho)

if isvector(A)
	p = length(A);
	p1 = p-1;
	A1 = [A(:)'; eye(p1) zeros(p1,1)]; % VAR coefficients for 1-lag problem
else
	[n,n1,p] = size(A);
	assert(n1 == n,'VAR/VMA coefficients matrix has bad shape');
	pn1 = (p-1)*n;
	A1 = [reshape(A,n,p*n); eye(pn1) zeros(pn1,n)]; % VAR coefficients for 1-lag problem
end

% calculate spectral norm

rho = max(abs(eig(A1,'nobalance'))); % v2.0 - don't balance!

if nargin < 2 || isempty(newrho)
    assert(nargout <= 1,'too many output parameters');
    out1 = rho;                     % spectral norm
else
    out1 = var_decay(A,newrho/rho); % adjusted coefficients
    out2 = rho;                     % previous value of spectral norm
end
