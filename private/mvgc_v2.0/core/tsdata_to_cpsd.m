function [S,f,fres] = tsdata_to_cpsd(X,mtaper,fs,nwin,sparms,fres,autospec,verb)

afreq = nargin < 3 || isempty(fs);

if afreq, fs = 2*pi; end % default to angular frequency

if nargin < 4, nwin   = []; end
if nargin < 5, sparms = []; end

[nvars,nobs,ntrials] = size(X);

% Note that for multi-trial data we are already averaging over trials,
% so we assume that longer windows - possibly the entire length of the trial
% (with no overlap) - are appropriate.

if isempty(nwin)
	if ntrials > 1
		nwin = nobs;          % sensible default
	else
		nwin = floor(nobs/8); % sensible default
	end
else
	if nwin < 0 % negative: -nwin is a *factor* of trial length (-1 for whole trial)
		nwin = floor(-nwin*nobs);
	end
end
assert(nwin > 0 && nwin <= nobs && isint(nwin),'Window size must be a positive integer not larger than the number of observations');
wins = floor(nobs/nwin);

if mtaper
	if isempty(sparms)
		tapers = [3 5]; % sensible default
	else
		tapers = sparms;
	end
else
	if wins > 1
		if isempty(sparms)
			noverlap = floor(nwin/2); % sensible default
		else
			if sparms < 0 % negative: -sparms is a *factor* of window length
				noverlap = floor(-sparms*nwin);
			else
				noverlap = sparms;
			end
		end
		assert(noverlap < nwin,'Overlap must be smaller than window');
	else
		noverlap = []; % N/A
	end
end

if nargin < 6 || isempty(fres)
	nfft = 2^(nextpow2(nwin));
	fres = nfft/2;
elseif fres < 0
	nfft = 2^(nextpow2(nwin)-fres); % -fres is padding power!
	fres = nfft/2;
else
	nfft = 2*fres;
end

if nargin < 7 || isempty(autospec), autospec = false; end

if nargin < 8 || isempty(verb), verb = 1; end

if verb > 0
	if mtaper,   mstr = 'multi-taper'; else, mstr  = 'Welch';             end
	if autospec, astr = 'auto';        else, astr  = 'cross';             end
	if afreq,    fstr = 'angular';     else, fstr  = sprintf('%g Hz',fs); end
	fprintf('\nMethod               : %s (%s-spectra)\n', mstr,astr);
	fprintf('Sampling frequency   : %s\n',   fstr);
	fprintf('Variables            : %d\n',   nvars);
	fprintf('Trials               : %d\n',   ntrials);
	fprintf('Trial observations   : %d',     nobs);
	if ~afreq, fprintf(' (%g secs)\n',nobs/fs); else, fprintf('\n'); end
	fprintf('Windows              : %d\n',   wins);
	fprintf('Window observations  : %d',   nwin);
	if ~afreq, fprintf(' (%g secs)\n',nwin/fs); else, fprintf('\n'); end
	if mtaper
		fprintf('Tapers               : [%d %d]\n', tapers(1),tapers(2));
	else
		if isempty(noverlap)
			if isempty(sparms)
				fprintf('Overlap              : N/A\n');
			else
				fprintf('Overlap              : ignored\n');
			end
		else
			fprintf('Overlap              : %d (',noverlap);
			if ~afreq, fprintf('%g secs, ',noverlap/fs); end
			fprintf('%g%% of window)\n',100*noverlap/nwin);
		end
	end
	fprintf('Frequency resolution : %d\n\n', fres);
end

if isempty(fs), fs = 2*pi; end % default to angular frequency

h = fres+1;
if autospec
	S = zeros(h,nvars);
else
	S = zeros(h,nvars,nvars);
end
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
			if autospec
				for i = 1:nvars
					S(:,i) = S(:,i) + 2*mean(conj(XFT(:,:,i)).*XFT(:,:,i),2);
				end
			else
				for i = 1:nvars
					for j = 1:nvars
						S(:,i,j) = S(:,i,j) + 2*mean(conj(XFT(:,:,i)).*XFT(:,:,j),2);
					end
				end
			end
		end
		if verb > 1, fprintf(' done\n'); end
	end

else

	X = permute(X,[2 1 3]);
	for r = 1:ntrials
		if verb > 1, fprintf('trial %d of %d ...',r,ntrials); end
		if autospec
			for i = 1:nvars
				S(:,i) = S(:,i) + pwelch(X(:,i,r),nwin,noverlap,nfft,fs);
			end
		else
			for i = 1:nvars
				S(:,:,i) = S(:,:,i) + cpsd(X(:,i,r),X(:,:,r),nwin,noverlap,nfft,fs);
			end
		end
		if verb > 1, fprintf(' done\n'); end
	end

end

if autospec
	S = (fs/2)*S/ntrials;
else
	S = (fs/2)*permute(S/ntrials,[3 2 1]);

	% NOTE: to extract the autospectra:
	%
	% n = nvars*nvars;
	% d = (1:nvars+1:n)+(0:fres)'*n; % all the diagonals!
	% Sas = S(d);

end

function XFT = mtfft(X,tapers,nfft,fs)

[nobs, nx] = size(X);            % size of data
[m1,nt]    = size(tapers);       % size of tapers
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
