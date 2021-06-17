function [mosvc,r] = tsdata_to_sssvc(X,h,r,plotm)

% plotm = []      - don't plot
% plotm = n       - Matlab plot to figure n (if zero, use next)
% plotm = string  - Gnuplot terminal (may be empty)

[n,m,N] = size(X);

assert(all(isint(h(:))),'past/future horizon must be a 2-vector or a scalar positive integer');
if isscalar(h)
    p = h;    f = h;
elseif isvector(h) && length(h) == 2
    p = h(1); f = h(2);
else
    error('past/future horizon must be a 2-vector or a scalar positive integer');
end
assert(p+f < m,'past/future horizon too large (or not enough data)');
rmax = n*min(p,f);

if nargin < 3 || isempty(r), r  = rmax;  end

assert(isscalar(r) && isint(r) && r >= 0 && r <= rmax,'model order must be empty, or a positive integer <= n*min(p,f) = %d',rmax);

X = demean(X); % no constant term (don't normalise!)

mp  = m-p;
mp1 = mp+1;
mf  = m-f;
mh  = mp1-f; % m-p-f+1

M  = N*mp;
M1 = N*mp1;
Mh = N*mh;

Xf = zeros(n,f,mh,N);
for k = 1:f
    Xf(:,k,:,:) = X(:,p+k:mf+k,:);
end
Xf = reshape(Xf,n*f,Mh);

XP = zeros(n,p,mp1,N);
for k = 0:p-1
    XP(:,k+1,:,:) = X(:,p-k:m-k,:);
end
Xp = reshape(XP(:,:,1:mh,:),n*p,Mh);
XP = reshape(XP,n*p,M1);

[Wf,cholp] = chol((Xf*Xf')/Mh,'lower');
assert(cholp == 0,'forward weight matrix not positive definite');

[Wp,cholp] = chol((Xp*Xp')/Mh,'lower');
assert(cholp == 0,'backward weight matrix not positive definite');

BETA = Xf/Xp; % 'OH' estimate: regress future on past
assert(all(isfinite(BETA(:))),'subspace regression failed');

[~,S] = svd(Wf\BETA*Wp); % SVD of CCA-weighted OH estimate

sval = diag(S);    % the singular values
df   = 2*n*(1:r)'; % number of free parameters (Hannan & Deistler, see also Bauer 2001) ... or r*r+2*n*r ???
svc  = -log(1-[sval(2:end);0]) + df*(log(Mh)/Mh); % Bauer's Singular Value Criterion

morder = (0:r)';
[~,idx] = min(svc); mosvc = morder(idx);

if ~isempty(plotm) % we're going to plot

	mo = (1:r)';

	if mosvc == r, wsvc = '*'; else wsvc = ''; end

	gap = 0.05;
	ssvc = gap+(1-gap)*(svc-min(svc))/(max(svc)-min(svc));

	if ischar(plotm) % Gnuplot

		gpstem = fullfile(tempdir,'sssvc');
		gpdat = [mo ssvc sval];
		gp_write(gpstem,gpdat);

		gp = gp_open(gpstem,plotm,[Inf,0.6]);

		fprintf(gp,'datfile = "%s.dat"\n',gpstem);

		fprintf(gp,'\nset grid\n');
		fprintf(gp,'set xr[0:%g]\n',r);
		fprintf(gp,'set xlabel "SS dimension"\n');

		fprintf(gp,'\nset multiplot title "SS SVC model order selection (CCA, max = %d)\\\n" layout 2,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n',r);

		fprintf(gp,'\nset title "Singular value criterion (SVC)"\n');
		fprintf(gp,'set ytics 0.2\n');
		fprintf(gp,'set yr[0:1.05]\n');
		fprintf(gp,'set ylabel "SVC (scaled)"\n');
		fprintf(gp,'set key top right Left rev\n');
		fprintf(gp,'plot \\\n');
		fprintf(gp,'datfile u 1:2 w linespoints pt 6 ps 1.4 t "SVC (opt = %2d%c)"\n',mosvc,wsvc);

		fprintf(gp,'\nset title "Singular values"\n');
		fprintf(gp,'unset logs y\n');
		fprintf(gp,'set ytics auto format ''%% h''\n');
		fprintf(gp,'set yr [0:*]\n');
		fprintf(gp,'set ytics 0.2 nomirror\n');
		fprintf(gp,'set ylabel "singular value"\n');
		fprintf(gp,'plot datfile u 1:3 w boxes fs solid 0.25 not\n');

		fprintf(gp,'\nunset multiplot\n');

		gp_close(gp,gpstem,plotm);

	else % Matlab

		if plotm == 0, figure; else, figure(plotm); end; clf;

		xlims = [0 r];

		subplot(2,1,1);
		plot(mo,ssvc,'o-');
		grid on
		title('Singular value criterion (SVC)');
		ylabel('SVC (scaled)');
		xlabel('SS dimension');
		legend(sprintf('SVC (opt = %d%c)',mosvc,wsvc));
		xlim(xlims);
		ylim([0 1+gap]);

		subplot(2,1,2); % SVC
		bar(mo,sval,1.01,'FaceColor',[0.65 0.75 1]);
		xlim(xlims);
		xlabel('SS dimension');
		ylabel('singular value');
		title('singular values');

		axes('Units','Normal');
		h = title(sprintf('SS SVC model order selection (CCA, max = %d)\n\n',r),'FontSize',13);
		set(gca,'visible','off')
		set(h,'visible','on')
	end
end
