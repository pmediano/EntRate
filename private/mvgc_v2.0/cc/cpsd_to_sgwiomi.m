function C = cpsd_to_sgwiomi(S,group)

[n,n1,h] = size(S);
assert(n1 == n,'CPSD must be square at each frequency');

g = check_group(group,n);

LDS = zeros(h,1);
for k = 1:h
    LDS(k)  = logdet(S(:,:,k));
end

C = nan(g,h);
for a = 1:g
	ga = group{a};
	gb = 1:n; gb(ga) = [];
    for k = 1:h
		C(a,k) = logdet(S(ga,ga,k)) + logdet(S(gb,gb,k)) - LDS(k);
	end
end
