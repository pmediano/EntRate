function L = ss_likelihood(y,A,C,Q,R,S)

[n,T]   = size(y);
[r,r1]  = size(A); assert(r1 == r);
[n1,r1] = size(C); assert(n1 == n && r1 == r);
[n1,n2] = size(R); assert(n1 == n && n2 == n);

assert(maxabs(R-R') < eps);
R = symmetrise(R);
assert(isposdef(R));

if nargin < 6 || isempty(S); % innovations form: Q is K, R is SIG
    [r1,n1] = size(Q); assert(r1 == r && n1 == n)
    S = Q*R;
    E = Q*chol(R,'lower');
    Q = E*E';
else
    [r1,r2] = size(Q); assert(r1 == r && r2 == r);
    [r1,n1] = size(S); assert(r1 == r && n1 == n);
end
assert(isequal(Q,Q'));
Q = symmetrise(Q);

try
    P = dlyap(A,Q); % P_1 = A*P_1*A' + Q (assuming stationarity)
catch except
    error(except.message);
end
assert(maxabs(P-P') < eps);
P = symmetrise(P);
assert(isposdef(P));

z = zeros(r,1); % z_t = E[x_{t|t-1}] : z_1 = 0

% "naive" Kalman filter (square root version would be better)
%
% V_t = E[e_t*e_t']
% K_t = Kalman gain
% e_t = y_t - C*z_t
% P_t = E[(x_t - x_{t|t-1})*(x_t - x_{t|t-1})']

L = 0; % log-likelihood
for t = 1:T
    V = C*P*C' + R;              % should be posdef
    I = V\eye(n);                % information matrix
    K = (A*P*C'+S)*I;            % gain matrix
    e = y(:,t) - C*z;            % prediction error

    L = L + (logdet(V) + e'*I*e)/T; % the log-likelihood

    P = A*P*A' + Q - K*V*K';     % theoretically positive-definite, numerically not so much...
    P = (P+P')/2;                % symmetrise for good measure (shouldn't be necessary)
    z = A*z + K*e;               % next state prediction
end

% Can force negative eigenvalues > 0, ensuring posdef... but is it a good idea?
% Does increase computation substantially...
%
% [U,D] = eig(P);
% D(D < 0) = eps;
% P = U*D*U';
