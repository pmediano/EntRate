%% plot_spw
%
% Plot spectral pairwise quantities on a grid
%
% <matlab:open('plot_spw.m') code>
%
%% Syntax
%
%     plot_spw(P,lam)
%
%% Arguments
%
% _input_
%
%     P          matrix (or cell vector of matrices) of spectral pairwise quantities
%     lam        vector of frequencies or sampling rate
%     ptitle     plot title (optional)
%
%% Description
%
% Plot pairwise spectral quantities in |P|, a 3-dim numerical matrix with
% first index representing target ("to"), second index source ("from")
% quantities and third index frequencies - typically spectral causalities (see
% e.g. <autocov_to_spwcgc.html |autocov_to_spwcgc|>).
%
%% See also
%
% <autocov_to_spwcgc.html |autocov_to_spwcgc|> |
% <mvgc_demo.html |mvgc_demo|> |
% <sfreqs.html |sfreqs|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_spw(P,lam,ptitle,plotm)

% plotm = n       - Matlab plot to figure n (if empty or zero, use next)
% plotm = string  - Gnuplot to terminal in string (may be empty)

if nargin < 4 || isempty(plotm), plotm = 0; end

if iscell(P)
	assert(isvector(P),'spectral quantities must be a (cell vector of) 3-dim matrix with the first two dims square');
else
	PP = P;
	clear P;
	P{1} = PP;
end
[n,~,h] = size(P{1});
S = length(P);

Pm = zeros(S,1);
Px = zeros(S,1);
for s = 1:S
	assert(ndims(P{s}) == 3 && size(P{s},2) == n,'spectral quantities must be 3-dim matrices with the first two dims square');
	if ~isreal(P{s}), fprintf(2,'WARNING(plot_spw): there were complex values!\n'); end
	P{s}(imag(P{s}) ~= 0) = NaN; % v2.0 - deal with complex values
	Pm(s) = min(P{s}(:));
	Px(s) = max(P{s}(:));
end

fres = h-1;

if nargin < 2 || isempty(lam)
	lam = sfreqs(fres)'; % normalised to [0,pi]
	xlab = 'normalised frequency';
else
	assert(isvector(lam),'frequencies must be supplied as a vector, or else a scalar sample frequency');
	if isscalar(lam)
		lam = sfreqs(fres,lam)'; % assume lam is the sample rate
	else
		assert(length(lam) == h,'frequency vector must match spectral quantities');
	end
	xlab = 'frequency (Hz)'; % assumed
end

xlims = [lam(1) lam(end)];
ylims = [min(Pm) 1.1*max(Px)];

if ischar(plotm) % Gnuplot

	fprintf(2,'WARNING: plot_spw: Gnuplot plotting not yet implemented\n');

else

	if plotm == 0, figure; else, figure(plotm); end; clf;

	dco = get(0,'DefaultAxesColorOrder'); % save default colour order
	set(0, 'DefaultAxesColorOrder',[0 0 0.8; 0.8 0 0; 0 0.8 0; dco]); % set to blue, red, green

	Pij = zeros(h,S);
	k = 0;
	for i = 1:n
		for j = 1:n
			k = k+1;
			if i ~= j
				subplot(n,n,k);
				for s = 1:S
					Pij(:,s) = squeeze(P{s}(i,j,:));
				end
				h = plot(lam,Pij);
				axis('square');
				xlim(xlims);
				ylim(ylims);
				xlabel(xlab);
				ylabel(sprintf('%d -> %d',j,i));
			end
		end
	end
	set(0, 'DefaultAxesColorOrder',dco); % restore default colour order

	if nargin > 2 && ~isempty(ptitle);
		axes('Units','Normal');
		h = title(sprintf('%s\n\n',ptitle));
		set(gca,'visible','off')
		set(h,'visible','on')
	end
end
