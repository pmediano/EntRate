function C = cpsd_to_spwcmi(S)

[n,n1,h] = size(S);
assert(n1 == n,'CPSD must be square at each frequency');

LDS = zeros(h,1);
for k = 1:h
    LDS(k)  = logdet(S(:,:,k));
end

LDSI = zeros(n,h);
for i = 1:n
    oi = 1:n; oi(i) = [];
    for k = 1:h
        LDSI(i,k) = logdet(S(oi,oi,k));
    end
end

C = nan(n,n,h);

for i = 1:n
	for j = i+1:n
        oij = 1:n; oij([i j]) = [];
        for k = 1:h
            C(i,j,k) = LDSI(i,k) + LDSI(j,k) - logdet(S(oij,oij,k)) - LDS(k);
        end
        C(j,i,:) = C(i,j,:);
    end
end
