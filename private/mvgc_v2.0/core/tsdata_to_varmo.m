function [moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X,p,regmode,alpha,pacf,plotm,verb)

% plotm = []      - don't plot
% plotm = n       - Matlab plot to figure n (if zero, use next)
% plotm = string  - Gnuplot terminal (may be empty)

if nargin < 4 || isempty(alpha), alpha = [0.01 0.05]; end
if nargin < 5 || isempty(pacf),  pacf  = true;        end
if nargin < 6,                   plotm = [];          end
if nargin < 7 || isempty(verb),  verb  = 0;           end

if isempty(plotm), pacf = false; end % no point if we're not going to display

if isscalar(alpha), alpha = [alpha alpha]; end

n = size(X,1);

% get log-likelihoods, pacfs at VAR model orders 1 .. p

[LL,Lk,Lm,R,RM] = tsdata_to_varll(X,p,regmode,pacf,verb>1); % log-likelihood, #free parameters. effective #observations

% calculate information criteria

[aic,bic,hqc] = infocrit(LL,Lk,Lm); % Akaike, Schwarz' Bayesian, Hannan-Quinn

% calculate optimal model orders according to information criteria (note: NaNs are ignored)

morder = (0:p)';
[~,idx] = min(aic); moaic = morder(idx);
[~,idx] = min(bic); mobic = morder(idx);
[~,idx] = min(hqc); mohqc = morder(idx);

% sequential likelihood ratio F-test for optimal model order (Lutkephol)

lambda = 2*diff(Lm.*LL); % LR test statistics
df = n*n;            % degrees of freedom (number of free parameters)
df2 = Lm(2:end)-n*(1:p)'-1;
lrpval = 1-fcdf(lambda/df,df,df2);
lrpval(df2 <= 0) = Inf; % hack (cough) - give up if F-test gives up
lralpha = alpha(1)/p;   % Bonferroni adjustment to significance level
hit = false;
for k = p:-1:1
    if lrpval(k) < lralpha, hit = true; break; end
end
if hit, molrt = k; else, molrt = 0; end

if verb > 0
    fprintf('\nBest model orders\n');
    fprintf('-----------------\n\n');
    fprintf('AIC : %2d',moaic); if moaic == p, fprintf(' *'); end; fprintf('\n');
    fprintf('BIC : %2d',mobic); if mobic == p, fprintf(' *'); end; fprintf('\n');
    fprintf('HQC : %2d',mohqc); if mohqc == p, fprintf(' *'); end; fprintf('\n');
    fprintf('LRT : %2d',molrt); if molrt == p, fprintf(' *'); end; fprintf('\n\n');
end

if ~isempty(plotm) % we're going to plot

	mo = morder;

	if moaic == p, waic = '*'; else waic = ''; end
	if mobic == p, wbic = '*'; else wbic = ''; end
	if mohqc == p, whqc = '*'; else whqc = ''; end
	if molrt == p, wlrt = '*'; else wlrt = ''; end

	gap = 0.05;
	saic = gap+(1-gap)*(aic-min(aic))/(max(aic)-min(aic));
	sbic = gap+(1-gap)*(bic-min(bic))/(max(bic)-min(bic));
	shqc = gap+(1-gap)*(hqc-min(hqc))/(max(hqc)-min(hqc));

	lmin = eps; lpval = lrpval; lpval(lpval<=0) = lmin;

	if pacf
		ccalpha = alpha(2)/(n*n*p);             % Bonferroni correction on significance levels
		Rcrit = norminv(1-ccalpha/2)./sqrt(RM); % 2-sided test
		R(:,:,1) = nan(n);                      % lag-zero are trivial, don't display
		Rlim = 1.1*nanmax(abs(R(:)));           % for display
	end

	if ischar(plotm) % Gnuplot

		gpstem = fullfile(tempdir,'varmo');

		gpdat = [mo saic sbic shqc [NaN;lpval]];

		if pacf
			RR = zeros(p+1,n*n);
			k = 1;
			for i = 1:n
				for j = 1:n
					RR(:,k) = R(i,j,:);
					k = k+1;
				end
			end

			gpdat = [gpdat RR Rcrit];
		end

		gp_write(gpstem,gpdat);

		gp = gp_open(gpstem,plotm,[Inf,0.6]);

		fprintf(gp,'datfile = "%s.dat"\n',gpstem);

		fprintf(gp,'\nset grid\n');
		fprintf(gp,'set xr[0:%g]\n',p);
		fprintf(gp,'set xlabel "AR lags"\n');

		if pacf
			fprintf(gp,'\nset multiplot title "VAR model order selection (%s, max = %d)\\\n" layout 3,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n',regmode,p);
		else
			fprintf(gp,'\nset multiplot title "VAR model order selection (%s, max = %d)\\\n" layout 2,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n',regmode,p);
		end

		fprintf(gp,'\nset title "Information criteria"\n');
		fprintf(gp,'set ytics 0.2\n');
		fprintf(gp,'set yr[0:1.05]\n');
		fprintf(gp,'set ylabel "information criterion (scaled)"\n');
		fprintf(gp,'set key top right Left rev\n');
		fprintf(gp,'plot \\\n');
		fprintf(gp,'datfile u 1:2 w linespoints pt 6 ps 1.4 t "AIC (opt = %2d%c)", \\\n',moaic,waic);
		fprintf(gp,'datfile u 1:3 w linespoints pt 6 ps 1.4 t "BIC (opt = %2d%c)", \\\n',mobic,wbic);
		fprintf(gp,'datfile u 1:4 w linespoints pt 6 ps 1.4 t "HQC (opt = %2d%c)"\n',mohqc,whqc);

		fprintf(gp,'\nset title "Likelihood ratio F-test (alpha = %.2g)"\n',alpha(1));
		fprintf(gp,'set ytics auto format ''%%.0e'' nomirror\n');
		fprintf(gp,'set yr[*:100]\n');
		fprintf(gp,'set ylabel "p-value"\n');
		fprintf(gp,'set logs y\n');
		fprintf(gp,'#set key top left Left rev\n');
		fprintf(gp,'set arrow 1 from graph 0, first %g to graph 1,first %g lt 2 nohead\n',lralpha,lralpha);
		fprintf(gp,'plot datfile u 1:5 w linespoints pt 6 ps 1.4 t "LRT (opt = %2d%c)"\n',molrt,wlrt);

		if pacf
			fprintf(gp,'\nset title "Partial autocorrelation (alpha = %.2g)"\n',alpha(2));
			fprintf(gp,'set format\n');
			fprintf(gp,'unset logs y\n');
			fprintf(gp,'unset key\n');
			fprintf(gp,'unset arrow 1\n');
			fprintf(gp,'set ytics auto mirror\n');
			fprintf(gp,'set yr[-%g:%g]\n',Rlim,Rlim);
			fprintf(gp,'set ylabel "correlation coefficient"\n');
			fprintf(gp,'plot \\\n');
			for k = 1:n*n
				fprintf(gp,'datfile u 1:%d w points pt 6 ps 1.2 not, \\\n',5+k);
			end
			fprintf(gp,'datfile u 1:(+($%d)) w lines lt 2 not, \\\n',6+n*n);
			fprintf(gp,'datfile u 1:(-($%d)) w lines lt 2 not\n',6+n*n);
		end

		fprintf(gp,'\nunset multiplot\n');

		gp_close(gp,gpstem,plotm);

	else % Matlab

		if plotm == 0, figure; else, figure(plotm); end; clf;

		xlims = [0 p];

		if pacf, subplot(3,1,1); else, subplot(2,1,1); end
		plot(mo,[saic sbic shqc],'o-');
		grid on
		title('Information criteria');
		ylabel('information criterion (scaled)');
		xlabel('AR lags');
		legend(sprintf('AIC (opt = %d%c)',moaic,waic),sprintf('BIC (opt = %d%c)',mobic,wbic),sprintf('HQC (opt = %d%c)',mohqc,whqc));
		xlim(xlims);
		ylim([0 1+gap]);

		if pacf, subplot(3,1,2); else, subplot(2,1,2); end
		semilogy(mo,[NaN;lpval],'o-');
		yline(lralpha,'r');
		grid on
		title(sprintf('Likelihood ratio F-test (\\alpha = %g)',alpha(1)));
		ylabel('p-value');
		xlabel('AR lags');
		legend(sprintf('LRT (opt = %d%c)',molrt,wlrt),'location','southeast');
		xlim(xlims);
		ylim([lmin 100]);

		if pacf
			subplot(3,1,3);
			hold on
			for i = 1:n
				for j = 1:n
					plot(mo,squeeze(R(i,j,:)),'o');
				end
			end
			plot(mo,[Rcrit -Rcrit],'r');
			hold off
			grid on
			title(sprintf('Partial autocorrelation (\\alpha = %g)',alpha(2)));
			ylabel('correlation coefficient');
			xlabel('AR lags');
			xlim(xlims);
			ylim([-Rlim Rlim]);
		end

		axes('Units','Normal');
		h = title(sprintf('VAR model order selection (%s, max = %d)\n\n',regmode,p),'FontSize',13);
		set(gca,'visible','off')
		set(h,'visible','on')
	end
end
