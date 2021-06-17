function [S,f,fres] = xspectrum(X,mtaper,fs,winfac,tapover,fres,verb)

afreq = nargin < 3 || isempty(fs);

if afreq, fs = 2*pi; end % default to angular frequency

if nargin < 4 || isempty(winfac), winfac = 8; end % sensible default

if nargin < 5 || isempty(tapover)
	if mtaper, tapers = [3 5];   else, noverlap = [];      end % sensible defaults
else
	if mtaper, tapers = tapover; else, noverlap = tapover; end
end

[nvars,nobs,ntrials] = size(X);

nwin = round(nobs/winfac);

if nargin < 6 || isempty(fres)
	nfft = 2^(nextpow2(nwin));
	fres = nfft/2;
elseif fres <= 0
	nfft = 2^(nextpow2(nwin)-fres); % -fres is padding power!
	fres = nfft/2;
else
	nfft = 2*fres;
end

assert(nwin <= nobs,'Window too large');
if nwin > nfft
	fprintf(2,'WARNING: resolution too low or window too large\n');
end

if nargin < 7 || isempty(verb), verb = 0; end

if verb > 0
	if mtaper, mstr = 'multi-taper'; else, mstr = 'Welch'; end
	if afreq, fstr  = 'angular'; else, fstr  = sprintf('%g Hz',fs);	end
	fprintf('\nMethod               : %s\n', mstr);
	fprintf('Sampling frequency   : %s\n',   fstr);
	fprintf('Variables            : %d\n',   nvars);
	fprintf('Trials               : %d\n',   ntrials);
	if ~afreq, fprintf('Trial time           : %g secs\n',nobs/fs); end
	fprintf('Trial observations   : %d\n',   nobs);
	if ~afreq, fprintf('Window time          : %g secs\n',nwin/fs); end
	fprintf('Window observations  : %d\n',   nwin);
	fprintf('Frequency resolution : %d\n\n', fres);
end

if isempty(fs), fs = 2*pi; end % default to angular frequency

h = fres+1;
S = zeros(h,nvars,nvars);
f = linspace(0,fs/2,h)';

if mtaper

	sz = size(tapers);
	if sz(1) == 1 && sz(2) == 2
		tapers = dpss(nwin,tapers(1),tapers(2))*sqrt(fs);
	else
		assert(sz(1) == nwin,'seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
	end

	nwins  = floor(nobs/nwin);

	for r = 1:ntrials
		if verb > 1, fprintf('trial %d of %d ...',r,ntrials); end
		for w = 1:nwins
			Xwr = X(:,1+(w-1)*nwin:w*nwin,r);
			XFT = mtfft(dtx(Xwr'),tapers,nfft,fs);
			XFT = XFT(1:h,:,:);
			for i = 1:nvars
				for j = 1:nvars
					S(:,i,j) = S(:,i,j) + mean(conj(XFT(:,:,i)).*XFT(:,:,j),2);
				end
			end
		end
		if verb > 1, fprintf(' done\n'); end
	end
	S = fs*permute(S/(nwins*ntrials),[3 2 1]);

else

	X = permute(X,[2 1 3]);
	for r = 1:ntrials
		if verb > 1, fprintf('trial %d of %d ...',r,ntrials); end
		for i = 1:nvars
			Si = cpsd(X(:,i,r),X(:,:,r),nwin,noverlap,nfft,fs);
			S(:,:,i) = S(:,:,i) + Si;
		end
		if verb > 1, fprintf(' done\n'); end
	end
	S = (fs/2)*permute(S/ntrials,[3 2 1]);

end

function XFT = mtfft(X,tapers,nfft,fs)

[nobs, nx] = size(X);      % size of data
[m1,nt] = size(tapers); % size of tapers
assert(m1 == nobs,'length of tapers is incompatible with length of data');

tapers = tapers(:,:,ones(1,nx)); % add channel indices to tapers
X      = X(:,:,ones(1,nt));      % add taper indices to data
X      = permute(X,[1 3 2]);     % reshape data to get dimensions to match those of tapers
XT     = X.*tapers;              % product of data with tapers
XFT    = fft(XT,nfft)/fs;        % fft of projected data

function y = dtx(x)

nobs = size(x,1);
a = zeros(nobs,2);
a(1:nobs,1) = (1:nobs)/nobs;
a(1:nobs,2) = 1;
y = x - a*(a\x); % Remove best fit
