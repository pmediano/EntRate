function [F,rep] = vou_to_gc(A,V,x,y,verb)

% Calculate time-domain conditional Granger causality rate for a vector
% Ornstein-Uhlenbeck (VOU) process. Source/target/conditioning variables may be
% multivariate. Uses a state-space method which involves solving an associated
% continuous-time algebraic Riccati equation (CARE).
%
% A     - VOU coefficients matrix
% V     - VOU Wiener process covariance matrix
% x     - multi-index of target variable
% y     - multi-index of source variable
% verb  - verbosity flag (default: true)
%
% F     - Granger causality rate from y to x, conditional on other variables
% rep   - CARE report: negative value indicates an error condition
%         rep >=  0 : Frobenius norm of relative residual error of CARE solution
%         rep == -1 : CARE has eigenvalues on or close to the imaginary axis (fatal)
%         rep == -2 : CARE has no stabilising solutions (fatal)
%         rep == -3 : CARE solution not positive-definite (non-fatal)
%         rep == -4 : VOU covariance matrix not positive-definite (non-fatal)
%
% NOTE 1: The total number of variables (dimensionality of the system) is
%         identified from the size of the A, V matrices. All variables indices
%         that do NOT appear in the target x and source y multi-indices (which
%         must not overlap) are taken as conditioning variables.
%
% NOTE 2: The VOU is stable iff max(real(eig(A))) < 0; however, this function
%         will generally return a reasonable-looking answer even in the unstable
%         case. Whether that answer is meaningful (as a transfer entropy rate)
%         requires clarification :-/ See ref. (3) below.
%
% NOTE 3: This function requires the Matlab Control System Toolbox for the "care" function.
%
% EXAMPLE:
%
%    n = 8;                         % dimensionality of system
%    lambda = -1;                   % max. real eigenvalue (negative for stable system)
%    A = randn(n);                  % random VOU coefficients matrix
%    A = A - (max(real(eig(A)))-lambda)*eye(n); % adjust so that max. real eigenvalue = lambda
%    L = randn(n,40); V = (L*L')/n; % random positive-definite VOU covariance matrix
%    x = 1:3;                       % target variable indices
%    y = 4:5;                       % source variable indices
%    F = vou_to_gc(A,V,x,y);        % calculate conditional Granger causality rate
%    fprintf('Granger causality rate = %g nats per unit time\n',F);
%
% REFERENCES:
%
% (1) L. Barnett and A. K. Seth (2015): Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication.
% (2) Barnett and A. K. Seth (2016): Detectability of Granger causality for subsampled continuous-time neurophysiological processes, J. Neurosci. Methods 275.
% (3) L. Barnett (2017): Granger causality rate for a vector Ornstein-Uhlenbeck process (working notes).
%
% (C) Lionel Barnett, May 2017

if nargin < 5 || isempty(verb), verb = true; end

[n, n1]  = size(A); assert(n1 == n, 'VOU coeffcicients matrix must be square');
[n1,n2]  = size(V); assert(n1 == n2,'VOU covariance matrix must be square');
                    assert(n1 == n, 'VOU covariance matrix must be same size as coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),     'Some target indices out of range');
assert(all(y >=1 & y <= n),     'Some source indices out of range');
assert(isempty(intersect(x,y)), 'Source/target multi-indices overlap');

z = 1:n; z([x y]) = []; % indices of remaining (conditioning) variables
r = [x z];              % indices for the reduced system

% Solve the associated CARE (see refs. 1,3)

[P,~,~,rep] = care(A(y,y)',A(r,y)',V(y,y),V(r,r),V(r,y)');
if rep < 0 % show-stopper; otherwise, rep returns Frobenius norm of relative residual error of CARE solution
    if verb
        if rep == -1, fprintf(2,'ERROR: CARE has eigenvalues on or close to the imaginary axis\n'); end
        if rep == -2, fprintf(2,'ERROR: CARE has no stabilising solutions\n'); end
    end
    F = NaN;
    return
end

% Implement F = trace(V(x,x)\(A(x,y)*P*A(x,y)')), checking that covariance
% matrices are positive-semidefinite

[LP,cholp] = chol(P,'lower');
if cholp ~= 0
    rep = -3;
    if verb, fprintf(2,'WARNING: CARE solution not positive-definite\n'); end
    F = trace(V(x,x)\(A(x,y)*P*A(x,y)')); % fall-back
    return
end

[LV,cholp] = chol(V(x,x),'lower');
if cholp ~= 0
    rep = -4;
    if verb, fprintf(2,'WARNING: covariance matrix not positive-definite\n'); end
    F = trace(V(x,x)\(A(x,y)*P*A(x,y)')); % fall back
    return
end

W = LV\(A(x,y)*LP);

F = trace(W*W'); % Granger causality rate
