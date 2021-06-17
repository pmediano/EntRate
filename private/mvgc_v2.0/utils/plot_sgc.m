function plot_sgc(P,lam,ptitle,plotm)

% plotm = n       - Matlab plot to figure n (if empty or zero, use next)
% plotm = string  - Gnuplot terminal (may be empty)

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

	gpstem = fullfile(tempdir,'plot_sgc');

	PP = cell(n*(n-1),1);
	k = 0;
	for i = n:-1:1
		for j = 1:n
			if i ~= j
				k = k+1;
				PP{k}(:,1) = lam;
				for s = 1:S
					PP{k}(:,s+1) = squeeze(P{s}(i,j,:));
				end
			end
		end
	end

	gp_write(gpstem,PP);

	gp = gp_open(gpstem,plotm,[Inf Inf]);

	fprintf(gp,'datfile = "%s.dat"\n\n',gpstem);

	fprintf(gp,'set style line 1 lt 1 lc rgb "#0000AA" lw 2\n');
	fprintf(gp,'set style line 2 lt 1 lc rgb "#AA0000" lw 2\n');

	fprintf(gp,'\nset xlab "frequency (Hz)"\n');
	fprintf(gp,'set xr[%g:%g]\n',xlims(1),xlims(2));
	fprintf(gp,'set yr[%g:%g]\n',ylims(1),ylims(2));

	if nargin > 2 && ~isempty(ptitle)
		fprintf(gp,'\nset multiplot title "%s\\n\\n" layout %d,%d rowsfirst\n',ptitle,n,n);
	else
		fprintf(gp,'\nset multiplot layout %d,%d rowsfirst\n',n,n);
	end

	k = 0;
	for i = n:-1:1
		for j = 1:n
			if i == j
				fprintf(gp,'\nset multiplot next\n');
			else
				fprintf(gp,'\nset ylab "%d -> %d"\n',j,i);
				fprintf(gp,'plot datfile i %d u 1:2 with lines ls 1 not, datfile i %d u 1:3 with lines ls 2 not\n',k,k);
				k = k+1;
			end
		end
	end

	fprintf(gp,'\nunset multiplot\n');

	gp_close(gp,gpstem,plotm);

else

	if plotm == 0, figure; else, figure(plotm); end; clf;

	dco = get(0,'DefaultAxesColorOrder'); % save default colour order
	set(0, 'DefaultAxesColorOrder',[0 0 0.8; 0.8 0 0; 0 0.8 0; dco]); % set to blue, red, green

	Pij = zeros(h,S);
	k = 0;
	for i = n:-1:1
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

	if nargin > 2 && ~isempty(ptitle)
		axes('Units','Normal');
		h = title(sprintf('%s\n\n',ptitle));
		set(gca,'visible','off')
		set(h,'visible','on')
	end
end
