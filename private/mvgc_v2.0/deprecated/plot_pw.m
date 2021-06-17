%% plot_pw
%
% Plot pairwise quantities on a colourmapped grid
%
% <matlab:open('plot_pw.m') code>
%
%% Syntax
%
%     plot_pw(P,cm)
%
%% Arguments
%
% _input_
%
%     P          square matrix of pairwise quantities
%     ptitle     plot title (optional)
%     cm         colour map (default: something soothing)
%     maxP       maximum value of P
%
%% Description
%
% Plot pairwise quantities in |P|, a 2-dim square numerical matrix with
% first index representing target ("to") and second index source ("from")
% quantities, typically causalities, p-values, significances, etc. (see
% e.g. <autocov_to_pwcgc.html |autocov_to_pwcgc|>). Diagonal entries
% are ignored. A colormap |cm| may be supplied.
%
%% See also
%
% <autocov_to_pwcgc.html |autocov_to_pwcgc|> |
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_pw(P,ptitle,cm,maxP,plotm)

% plotm = n       - Matlab plot to figure n (if empty or zero, use next)
% plotm = string  - Gnuplot to terminal in string (may be empty)

if nargin < 3 || isempty(cm),    cm    = flipud(bone); end
if nargin < 4 || isempty(maxP),  maxP  = 0;            end
if nargin < 5 || isempty(plotm), plotm = 0;            end

if iscell(P)
	n = numel(P);
	[nr,nc] = size(P);
	P = P(:);
	ptitle = ptitle(:);
	assert(iscell(ptitle) && numel(ptitle) >= n,'Not enough plot titles');
	if iscell(cm)
		assert(iscell(cm) && numel(cm) >= n,'Not enough colour maps');
		cm = cm(:);
	else
		cmtmp = cm;
		cm = cell(n,1);
		for i = 1:n, cm{i} = cmtmp; end
	end
	if isscalar(maxP)
		maxP = maxP*ones(n,1);
	else
		assert(numel(maxP) >= n,'Not enough max values');
		maxP = maxP(:);
	end

	for i = 1:n
		if isempty(P{i}), continue; end
		m = size(P{i},1);
		assert(ismatrix(P{i}) && size(P{i},2) == m,'input must be a square 2D matrix');
		if ~isreal(P{i}), fprintf(2,'WARNING(plot_pw): there were complex values!\n'); end
		P{i}(imag(P{i}) ~= 0) = NaN; % v2.0 - deal with complex values
		if maxP(i) == 0, maxP(i) = max(P{i}(:)); end
	end

	if ischar(plotm) % Gnuplot
		plot_pw_gp(P,ptitle,cm,maxP,nr,nc,plotm);
	else
		if plotm == 0, figure; else, figure(plotm); end; clf;
		for i = 1:n, plot_pw_mat(P{i},ptitle{i},cm{i},maxP(i),nr,nc,i); end
	end

else

	m = size(P,1);
	assert(ismatrix(P) && size(P,2) == m,'input must be a square 2D matrix');
	if ~isreal(P), fprintf(2,'WARNING(plot_pw): there were complex values!\n'); end
	P(imag(P) ~= 0) = NaN; % v2.0 - deal with complex values
	if maxP == 0, maxP = max(P(:)); end

	if ischar(plotm) % Gnuplot
		plot_pw_gp(P,ptitle,cm,maxP,nr,nc,plotm);
	else
		if plotm == 0, figure; else, figure(plotm); end; clf;
		plot_pw_mat(P,ptitle,cm,1,1,maxP);
	end

end

end

function plot_pw_gp(P,ptitle,cm,maxP,nr,nc,gpterm)

	n = numel(P);

	gpstem = fullfile(tempdir,'plotpw');

	gp_write_matrix(gpstem,P,false);

	gp = gp_open(gpstem,gpterm,1);

	fprintf(gp,'datfile = "%s.dat"\n',gpstem);
	fprintf(gp,'set palette defined ( 0 "#ffffff", 3 "#99bbcc", 10 "#000000" ) nops_allcF maxcolors 0\n'); % bone-like
	fprintf(gp,'#unset colorbox\n');
	fprintf(gp,'set xlab "from"\n');
	fprintf(gp,'set ylab "to"\n');

	fprintf(gp,'\nset multiplot layout %d,%d rowsfirst\n',nr,nc);

	for i = 1:n
		if isempty(P{i}), continue; end
		[r,c] = size(P{i});
		if r == c, fprintf(gp,'\nset size square\n'); end
		fprintf(gp,'set title "%s"\n',ptitle{i});
		fprintf(gp,'set xr[-0.5:%d-0.5]\n',c);
		fprintf(gp,'set yr[-0.5:%d-0.5]\n',r);
		fprintf(gp,'set xtics (%s)\n',gp_numtics(1:c,0:c-1));
		fprintf(gp,'set ytics (%s)\n',gp_numtics(r:-1:1,0:r-1));
		fprintf(gp,'set cbrange [0:%g]\n',maxP(i));
[i maxP(i)]
		fprintf(gp,'plot datfile i %d matrix with image not\n',i-1);
	end

	fprintf(gp,'\nunset multiplot\n');

	gp_close(gp,gpstem,gpterm);

end

function plot_pw_mat(P,ptitle,cm,maxP,nr,nc,i)

	if isempty(P), return; end

	if nargin > 4, subplot(nr,nc,i); end

	n = size(P,1);

	colormap(cm);
	imagesc(P,[0 maxP]);
	axis('square');
	xlabel('from');
	ylabel('to');
	set(gca,'XTick',1:n);
	set(gca,'XTickLabel',1:n);
	set(gca,'YTick',1:n);
	set(gca,'YTickLabel',1:n);
	title(ptitle)

end
