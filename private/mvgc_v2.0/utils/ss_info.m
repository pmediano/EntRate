function info = ss_info(A,C,K,V,report)

if nargin < 5 || isempty(report), report = 1; end % default: print out report

[r, r1] = size(A); assert(r1 == r);
[n, r1] = size(C); assert(r1 == r);
[r1,n1] = size(K); assert(n1 == n && r1 == r);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

info.error = uint32(0);

info.observ = n;
info.morder = r;

info.rhoA = max(abs(eig(A,'nobalance')));

B = A-K*C;
info.rhoB = max(abs(eig(B,'nobalance')));

info.acdec = ceil(0.5*log(eps)/log(max(info.rhoA,info.rhoB))); % so that autocov decays to < sqrt(eps), (probably ~ 1.5e-8)

if maxabs(triu(V,1)-triu(V',1)) > eps
    info.sigspd = 1; % not symmetric
else
    [~,cholp] = chol(V,'lower');
    if cholp > 0
        info.sigspd = 2; % symmetric, but not positive definite
    else
        info.sigspd = 0; % symmetric, positive definite
    end
end
info.mii = multiinfo(V);       % multi-information (generalised correlation)
if n > 1
  info.mmii = multiinfo(n,true); % multi-information for uniform random n x n correlation matrix, for comparison
else
  info.mmii = 0;
end

%{
if rand < 0.1 % test intermittent errors
    info.mii = -0.5;
    info.rhoB = 1;
end
%}

rhotol = sqrt(eps);

if     info.rhoA > 1+rhotol, info.error = bitset(info.error,1); % explosive
elseif info.rhoA > 1-rhotol, info.error = bitset(info.error,2); % unit root
end

if     info.rhoB > 1+rhotol, info.error = bitset(info.error,3); % explosive
elseif info.rhoB > 1-rhotol, info.error = bitset(info.error,4); % unit root
end

if     info.sigspd == 1,     info.error = bitset(info.error,5); % not symmetric
elseif info.sigspd == 2,     info.error = bitset(info.error,6); % not positive definite
end

if     info.mii < -1e-10,         info.error = bitset(info.error,7); % negative
end

if report == 1 % print out report

    fprintf('\nSS info:\n');

    fprintf('    observables       = %d\n',info.observ);

    fprintf('    model order       = %d\n',info.morder);

    fprintf('    AR spectral norm  = %.6f',info.rhoA);
    if      bitget(info.error,1), fprintf(2,'      ERROR: unstable (explosive)\n');
    elseif  bitget(info.error,2), fprintf(2,'      ERROR: unstable (unit root)\n');
    else    fprintf('      stable\n');
    end

    fprintf('    MA spectral norm  = %.6f',info.rhoB);
    if      bitget(info.error,3), fprintf(2,'      ERROR: not minimum phase (explosive)\n');
    elseif  bitget(info.error,4), fprintf(2,'      ERROR: not minimum phase (unit root)\n');
    else    fprintf('      minimum phase\n');
    end

    fprintf('    residuals covariance matrix');
    if     bitget(info.error,5), fprintf(2,'      ERROR: not symmetric\n');
    elseif bitget(info.error,6), fprintf(2,'      ERROR: not positive definite\n');
    else   fprintf('       symmetric, positive definite\n');
    end

    fprintf('    multi-information = %-.6f',info.mii);
    if     bitget(info.error,7), fprintf(2,'     ERROR: negative\n');
    else   fprintf('      uniform = %-.6f\n',info.mmii);
    end

    fprintf('    autocorr. decay   = %-7d\n',info.acdec);

    fprintf('\n');

elseif report > 1 % format error message(s) string

    if ~info.error, info.errmsg = ''; return; end % no errors to report

    info.nerrors = nnz(bitget(info.error,1:8)); % number of errors

    if info.nerrors > 1
        info.errmsg = 'SS ERRORS';
    else
        info.errmsg = 'SS ERROR';
    end

    if      bitget(info.error,1), info.errmsg = [info.errmsg sprintf(': AR spectral norm = %.6f - unstable (explosive)',info.rhoA)];
    elseif  bitget(info.error,2), info.errmsg = [info.errmsg sprintf(': AR spectral norm = %.6f - unstable (unit root)',info.rhoA)];
    end

    if      bitget(info.error,3), info.errmsg = [info.errmsg sprintf(': MA spectral norm = %.6f - not miniphase (explosive)',info.rhoB)];
    elseif  bitget(info.error,4), info.errmsg = [info.errmsg sprintf(': MA spectral norm = %.6f - not miniphase (unit root)',info.rhoB)];
    end

    if     bitget(info.error,5), info.errmsg = [info.errmsg ': res. cov. matrix not symmetric'];
    elseif bitget(info.error,6), info.errmsg = [info.errmsg ': res. cov. matrix not positive definite'];
    end

    if     bitget(info.error,7), info.errmsg = [info.errmsg sprintf(': multi-information = %.6f - negative',info.mii)];
    end

end
