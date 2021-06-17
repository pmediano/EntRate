function [moaic,mobic,mohqc,mosvc,molrt,ssmomax] = tsdata_to_ssmo(X,pf,morder,alpha,plotm,verb)

% plotm = []      - don't plot
% plotm = n       - Matlab plot to figure n (if zero, use next)
% plotm = string  - Gnuplot terminal (may be empty)

if nargin < 3,                   morder = [];    end
if nargin < 4 || isempty(alpha), alpha  = 0.01;  end
if nargin < 5,                   plotm = [];     end
if nargin < 6 || isempty(verb),  verb   = false; end

svconly = nargout < 3;

% get log-likelihoods

if svconly
	[~,~,~,~,svc,r] = tsdata_to_ssll(X,pf,morder,verb); % log-likelihood, #free parameters. effective #observations
	[~,idx] = min(svc);
	morder = (0:r)';
	moaic = morder(idx); % mosvc   !!!
	mobic = r;           % ssmomax !!!
	return
end

n = size(X,1);

[LL,Lk,Lm,sval,svc,r] = tsdata_to_ssll(X,pf,morder,verb); % log-likelihood, #free parameters. effective #observations
ssmomax = r;

% calculate information criteria

[aic,bic,hqc] = infocrit(LL,Lk,Lm); % Akaike, Schwarz' Bayesian, Hannan-Quinn

% calculate optimal model orders according to information criteria (note: NaNs are ignored)

morder = (0:r)';
[~,idx] = min(aic); moaic = morder(idx);
[~,idx] = min(bic); mobic = morder(idx);
[~,idx] = min(hqc); mohqc = morder(idx);
[~,idx] = min(svc); mosvc = morder(idx);

% sequential likelihood ratio F-test for optimal model order (Lutkephol)

lambda = [NaN;2*diff(Lm.*LL)]; % LR test statistics
rr = (1:r)';
df = 2*n;               % degrees of freedom (number of free parameters)
df2 = Lm-n*rr-1;
%size(df2)
%size(lambda)
lrpval = 1-fcdf(lambda./df,df,df2);
lrpval(df2 <= 0) = Inf; % hack (cough) - give up if F-test gives up
lralpha = alpha/r;      % Bonferroni adjustment to significance level
hit = false;
for k = r:-1:1
    if lrpval(k) < lralpha, hit = true; break; end
end
if hit, molrt = k; else, molrt = 0; end

if verb > 0
    fprintf('\nBest model orders\n');
    fprintf('-----------------\n\n');
    fprintf('AIC : %2d',moaic); if moaic == r, fprintf(' *'); end; fprintf('\n');
    fprintf('BIC : %2d',mobic); if mobic == r, fprintf(' *'); end; fprintf('\n');
    fprintf('HQC : %2d',mohqc); if mohqc == r, fprintf(' *'); end; fprintf('\n');
    fprintf('LRT : %2d',molrt); if molrt == r, fprintf(' *'); end; fprintf('\n');
    fprintf('SVC : %2d',mosvc); if mosvc == r, fprintf(' *'); end; fprintf('\n\n');
end

if ~isempty(plotm) % we're going to plot

	mo = (1:r)';

	if moaic == r, waic = '*'; else waic = ''; end
	if mobic == r, wbic = '*'; else wbic = ''; end
	if mohqc == r, whqc = '*'; else whqc = ''; end
	if molrt == r, wlrt = '*'; else wlrt = ''; end
	if mosvc == r, wsvc = '*'; else wsvc = ''; end

	gap = 0.05;
	saic = gap+(1-gap)*(aic-min(aic))/(max(aic)-min(aic));
	sbic = gap+(1-gap)*(bic-min(bic))/(max(bic)-min(bic));
	shqc = gap+(1-gap)*(hqc-min(hqc))/(max(hqc)-min(hqc));
	ssvc = gap+(1-gap)*(svc-min(svc))/(max(svc)-min(svc));

	lmin = eps; lpval = lrpval; lpval(lpval<=0) = lmin;

	if ischar(plotm) % Gnuplot

		gpstem = fullfile(tempdir,'ssmo');
		gpdat = [mo saic sbic shqc ssvc lpval sval];
		gp_write(gpstem,gpdat);

		gp = gp_open(gpstem,plotm,[Inf,0.6]);

		fprintf(gp,'datfile = "%s.dat"\n',gpstem);

		fprintf(gp,'\nset grid\n');
		fprintf(gp,'set xr[0:%g]\n',r);
		fprintf(gp,'set xlabel "SS dimension"\n');

		fprintf(gp,'\nset multiplot title "SS model order selection (CCA, max = %d)\\\n" layout 3,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n',ssmomax);

		fprintf(gp,'\nset title "Information criteria"\n');
		fprintf(gp,'set ytics 0.2\n');
		fprintf(gp,'set yr[0:1.05]\n');
		fprintf(gp,'set ylabel "information criterion (scaled)"\n');
		fprintf(gp,'set key top right Left rev\n');
		fprintf(gp,'plot \\\n');
		fprintf(gp,'datfile u 1:2 w linespoints pt 6 ps 1.4 t "AIC (opt = %2d%c)", \\\n',moaic,waic);
		fprintf(gp,'datfile u 1:3 w linespoints pt 6 ps 1.4 t "BIC (opt = %2d%c)", \\\n',mobic,wbic);
		fprintf(gp,'datfile u 1:4 w linespoints pt 6 ps 1.4 t "HQC (opt = %2d%c)", \\\n',mohqc,whqc);
		fprintf(gp,'datfile u 1:5 w linespoints pt 6 ps 1.4 t "SVC (opt = %2d%c)"\n',mosvc,wsvc);

		fprintf(gp,'\nset title "Likelihood ratio F-test (alpha = %.2g)"\n',alpha);
		fprintf(gp,'set ytics auto format ''%%.0e'' nomirror\n');
		fprintf(gp,'set yr[*:100]\n');
		fprintf(gp,'set ylabel "p-value"\n');
		fprintf(gp,'set logs y\n');
		fprintf(gp,'#set key top left Left rev\n');
		fprintf(gp,'set arrow 1 from graph 0, first %g to graph 1,first %g lt 2 nohead\n',lralpha,lralpha);
		fprintf(gp,'plot datfile u 1:6 w linespoints pt 6 ps 1.4 t "LRT (opt = %2d%c)"\n',molrt,wlrt);

		fprintf(gp,'\nset title "Singular values"\n');
		fprintf(gp,'unset logs y\n');
		fprintf(gp,'set ytics auto format ''%% h''\n');
		fprintf(gp,'set yr [0:*]\n');
		fprintf(gp,'set ytics 0.2 nomirror\n');
		fprintf(gp,'set ylabel "singular value"\n');
		fprintf(gp,'plot datfile u 1:7 w boxes fs solid 0.25 not\n');

		fprintf(gp,'\nunset multiplot\n');

		gp_close(gp,gpstem,plotm);

	else % Matlab

		if plotm == 0, figure; else, figure(plotm); end; clf;

		xlims = [0 r];

		subplot(3,1,1);
		plot(mo,[saic sbic shqc ssvc],'o-');
		grid on
		title('Information criteria');
		ylabel('information criterion (scaled)');
		xlabel('SS dimension');
		legend(sprintf('AIC (opt = %d%c)',moaic,waic),sprintf('BIC (opt = %d%c)',mobic,wbic),sprintf('HQC (opt = %d%c)',mohqc,whqc),sprintf('SVC (opt = %d%c)',mosvc,wsvc));
		xlim(xlims);
		ylim([0 1+gap]);

		subplot(3,1,2);
		semilogy(mo,lpval,'o-');
		yline(lralpha,'r');
		grid on
		title(sprintf('Likelihood ratio F-test (\\alpha = %g)',alpha));
		ylabel('p-value');
		xlabel('SS dimension');
		legend(sprintf('LRT (opt = %d%c)',molrt,wlrt),'location','southeast');
		xlim(xlims);
		ylim([lmin 100]);

		subplot(3,1,3); % SVC
		bar(mo,sval,1.01,'FaceColor',[0.65 0.75 1]);
		xlim(xlims);
		xlabel('SS dimension');
		ylabel('singular value');
		title('singular values');

		axes('Units','Normal');
		h = title(sprintf('SS model order selection (CCA, max = %d)\n\n',ssmomax),'FontSize',13);
		set(gca,'visible','off')
		set(h,'visible','on')
	end
end
