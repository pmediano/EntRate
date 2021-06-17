%% logdet
%
% Calculate logarithm of determinant for positive-definite matrices
%
% <matlab:open('logdet.m') code>
%
%% Syntax
%
%     L = logdet(V)
%
%% Arguments
%
% _input_
%
%     V          a positive-definite square, symmetric matrix
%
% _output_
%
%     L          logarithm of determinant of V
%
%% Description
%
% Returns the logarithm of determinant of the symmetric, positive-definite matrix |V| in |LD|.
% Essentially the same as |log(det(V))|, but avoids potential under/overflow, Does not check
% whether V is Hermitian, but returns |NaN| if V is not positive-definite.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function LD = logdet(V)

[L,cholp] = chol(V);
if cholp == 0
	LD = 2*sum(log(diag(L)));
else
	DV = det(V);
	if abs(imag(DV)) > sqrt(eps)
		LD = NaN;
	else % give it the benefit of the doubt...
		LD = log(real(DV));
	end
end
