function G = iss2cov(A,C,K,V)

L = chol(V,'lower');
KL = K*L;
U = C*chol(lyapslv('D',A,[],-KL*KL'),'lower');
G = chol(U*U'+L*L','lower');
G = G*G';
