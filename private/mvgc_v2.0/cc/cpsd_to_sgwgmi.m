function C = cpsd_to_sgwgmi(S,group)

[n,n1,h] = size(S);
assert(n1 == n,'CPSD must be square at each frequency');

[g,gsiz] = check_group(group,n);

LDS = zeros(h,1);
for k = 1:h
    LDS(k)  = logdet(S(:,:,k));
end

C = nan(g,h);

for a = 1:g
    goa = 1:n; goa(group{a}) = [];
    for k = 1:h
        C(a,k) = -LDS(k) - (gsiz(a)-1)*logdet(S(goa,goa,k));
    end
    for i = group{a}
        igoa = [i goa];
        for k = 1:h
            C(a,k) = C(a,k) + logdet(S(igoa,igoa,k));
        end
    end
end
