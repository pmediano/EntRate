function PBL = bandlimit(P,dim,fs,B)

% Assumption: P is defined from 0 - Nyqvist (inclusive) at regular spacing!

if nargin < 2 || isempty(dim), dim = 1; end

if nargin < 3 || isempty(fs)
    fs = 2*pi; % angular frequency
end

if nargin < 4, B = []; end

nd = ndims(P);
assert(dim >= 1 && dim <= nd,'bad dimension argument');

N = size(P,dim);

if isempty(B) % [0, Nyqvist]
    imin = 1;
    imax = N;
else
    bmsg = 'Frequency band must be an ascending 2-vector of frequencies in the range [0, Nyqvist]';
    assert(isnumeric(B) && isvector(B) && length(B) == 2,bmsg);
    fN = fs/2;
    fac = (N-1)/fN;
    if isnan(B(1))
        B(1) = 0;
        imin = 1;
    else
        assert(B(1) >= 0 && B(1) <= fN,bmsg);
        imin = 1+round(fac*B(1));
        if imin < 1, imin = 1; end
    end
    if isnan(B(2))
        B(2) = fN;
        imax = N;
    else
        assert(B(2) >= 0 && B(2) <= fN,bmsg);
        imax = 1+round(fac*B(2));
        if imax > N, imax = N; end
    end
    assert(B(1) < B(2),bmsg);
    assert(imin < imax,'frequency band limits too close together?');
end

od = 1:nd;    od(dim) = []; % other dimensions - all except dim
sz = size(P); sz(dim) = []; % other sizes      - all except dim
if isscalar(sz), sz = [sz 1]; end

P  = permute(P,[dim od]); % put quadrature dimension first
PBL = reshape(trapz(P(imin:imax,:))/(N-1),sz);
