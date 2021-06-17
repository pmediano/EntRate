function C = cpsd_to_sgwcmi(S,group)

[n,n1,h] = size(S);
assert(n1 == n,'CPSD must be square at each frequency');

g = check_group(group,n);;

LDS = zeros(h,1);
for k = 1:h
    LDS(k)  = logdet(S(:,:,k));
end

LDSG = zeros(g,h);
for a = 1:g
    goa = 1:n; goa(group{a}) = [];
    for k = 1:h
        LDSG(a,k) = logdet(S(goa,goa,k));
    end
end

C = nan(g,g,h);

for a = 1:g
	for b = a+1:g
        goab = 1:n; goab([group{a} group{b}]) = [];
        for k = 1:h
            C(a,b,k) = LDSG(a,k) + LDSG(b,k) - logdet(S(goab,goab,k)) - LDS(k);
        end
        C(b,a,:) = C(a,b,:);
    end
end
