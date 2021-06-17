function [moaic,mobic,mohqc,molrt] = tsdata_to_varmo_gp(X,p,regmode,alpha,pacf,plotit,verb,gpterm)

if nargin < 4 || isempty(alpha),   alpha   = 0.01;  end
if nargin < 5 || isempty(pacf),    pacf    = true;  end
if nargin < 6 || isempty(plotit),  plotit  = true;  end
if nargin < 7 || isempty(verb),    verb    = 0;     end
if nargin < 8 || isempty(gpterm),  gpterm  = '';    end

if ~plotit, pacf = false; end % no point if we're not going to display

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

lambda = 2*diff(LL); % LR test statistics
df = n*n;            % degrees of freedom (number of free parameters)
df2 = Lm(2:end)-n*(1:p)'-1;
lrpval = 1-fcdf(lambda/df,df,df2);
lrpval(df2 <= 0) = Inf; % hack (cough) - give up if F-test gives up
lralpha = alpha/p;   % Bonferroni correction on significance levels
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

if plotit

    mo = morder(2:end);
    xo = 0.25;

    if moaic == p, waic = '*'; else waic = ''; end
    if mobic == p, wbic = '*'; else wbic = ''; end
    if mohqc == p, whqc = '*'; else whqc = ''; end
    if molrt == p, wlrt = '*'; else wlrt = ''; end

    aic1 = aic(2:end);
    bic1 = bic(2:end);
    hqc1 = hqc(2:end);
    saic = 0.05+0.95*(aic1-min(aic1))/(max(aic1)-min(aic1));
    sbic = 0.05+0.95*(bic1-min(bic1))/(max(bic1)-min(bic1));
    shqc = 0.05+0.95*(hqc1-min(hqc1))/(max(hqc1)-min(hqc1));

    lmin = eps; lpval = lrpval; lpval(lpval<=0) = lmin;

    gpstem = fullfile(tempdir,'gptmp');
    gpdat = [mo saic sbic shqc lpval];
    if pacf

		R  = R(:,:,2:end);
		RM = RM(2:end);
		ccalpha = alpha/(n*n*p);                % Bonferroni correction on significance levels
		Rcrit = norminv(1-ccalpha/2)./sqrt(RM); % 2-sided test
        RR = zeros(p,n*n);
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

    gp = gp_open(gpstem,gpterm,[Inf,0.5]);

    fprintf(gp,'datfile = "%s.dat"\n',gpstem);

    fprintf(gp,'\nset grid\n');
    fprintf(gp,'set xr[%g:%g]\n',1-xo,p+xo);
    fprintf(gp,'set xlabel "lags"\n');

    if pacf
        fprintf(gp,'\nset multiplot layout 3,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n');
    else
        fprintf(gp,'\nset multiplot layout 2,1 margins 0.12,0.94,0.05,0.95 spacing 0.1\n');
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

	fprintf(gp,'\nset title "Likelihood ratio F-test (alpha = %.3f)"\n',alpha);
    fprintf(gp,'set ytics auto format ''%%.0e''\n');
    fprintf(gp,'set yr[*:100]\n');
    fprintf(gp,'set ylabel "p-value"\n');
    fprintf(gp,'set logs y\n');
    fprintf(gp,'#set key top left Left rev\n');
    fprintf(gp,'set arrow 1 from graph 0, first %g to graph 1,first %g lt 2 nohead\n',lralpha,lralpha);
    fprintf(gp,'plot datfile u 1:5 w linespoints pt 6 ps 1.4 t "LRT (opt = %2d%c)"\n',molrt,wlrt);

    if pacf
        fprintf(gp,'\nset title "Partial autocorrelation (alpha = %.3f)"\n',alpha);
        fprintf(gp,'set format\n');
        fprintf(gp,'unset logs y\n');
        fprintf(gp,'unset key\n');
        fprintf(gp,'unset arrow 1\n');
        fprintf(gp,'set ytics 0.5\n');
        fprintf(gp,'set yr[-1:1]\n');
        fprintf(gp,'set ylabel "correlation coefficient"\n');
        fprintf(gp,'plot \\\n');
        for k = 1:n*n
            fprintf(gp,'datfile u 1:%d w points pt 6 ps 1.2 not, \\\n',5+k);
        end
        fprintf(gp,'datfile u 1:(+($%d)) w lines lt 2 not, \\\n',6+n*n);
        fprintf(gp,'datfile u 1:(-($%d)) w lines lt 2 not\n',6+n*n);
    end

    fprintf(gp,'\nunset multiplot\n');

    gp_close(gp,gpstem,gpterm);

end
