function plot_gc(F,ptitle,cm,Fmax,plotm)

% plotm = n       - Matlab plot to figure n; if zero (default), use next Matlab figure
% plotm = string  - Gnuplot terminal (may be empty)

if nargin < 4 || isempty(Fmax),  Fmax  = 0; end
if nargin < 5 || isempty(plotm), plotm = 0; end
if nargin < 3 || isempty(cm)
	if ischar(plotm)
		cm = 'defined ( 0 "#ffffff", 3 "#99bbcc", 10 "#000000" ) nops_allcF maxcolors 0'; % bone-like
	else
		cm = flipud(bone);
	end
end

if iscell(F)
	assert(ismatrix(F),'GCs must be a cell matrix');
	[nr,nc] = size(F);
	assert(iscell(ptitle) && isequal(size(ptitle),[nr nc]),'Plot titles don''t match GCs');
	if isscalar(Fmax)
		Fmax = Fmax*ones(nr,nc);
	else
		assert(ismatrix(Fmax) && isequal(size(Fmax),[nr nc]),'Maximum values don''t match GCs');
	end

	for r = 1:nr
		for c = 1:nc
			if isempty(F{r,c}), continue; end
			assert(ismatrix(F{r,c}),'GCs must be matrices');
			if ~isreal(F{r,c}), fprintf(2,'WARNING: plot_gc: there were complex values!\n'); end
			F{r,c}(imag(F{r,c}) ~= 0) = NaN; % v2.0 - deal with complex values
			if Fmax(r,c) == 0, Fmax(r,c) = max(F{r,c}(:)); end
		end
	end

	if ischar(plotm) % Gnuplot
		plot_gc_gp(F,ptitle,cm,Fmax,plotm);
	else
		if plotm == 0, figure; else, figure(plotm); end; clf;
		i = 0;
		for r = 1:nr
			for c = 1:nc
				i = i+1;
				plot_gc_mat(F{r,c},ptitle{r,c},cm,Fmax(r,c),nr,nc,i);
			end
		end
	end

else

	assert(ismatrix(F),'GC must be a matrix');
	if ~isreal(F), fprintf(2,'WARNING: plot_gc: there were complex values!\n'); end
	F(imag(F) ~= 0) = NaN; % v2.0 - deal with complex values
	if Fmax == 0, Fmax = max(F(:)); end

	if ischar(plotm) % Gnuplot
		plot_gc_gp(F,ptitle,cm,Fmax,plotm);
	else
		if plotm == 0, figure; else, figure(plotm); end; clf;
		plot_gc_mat(F,ptitle,cm,Fmax);
	end

end

end

function plot_gc_gp(F,ptitle,cm,Fmax,gpterm)

	gpstem = fullfile(tempdir,'plot_gc');

	gp = gp_open(gpstem,gpterm);

	fprintf(gp,'set palette %s\n',cm);
	fprintf(gp,'unset colorbox\n');
	fprintf(gp,'set xlab "from"\n');
	fprintf(gp,'set ylab "to"\n');

	if iscell(F)
		[nr,nc] = size(F);
		fprintf(gp,'\nset multiplot layout %d,%d rowsfirst\n',nr,nc);
		for r = 1:nr
			for c = 1:nc
				if isempty(F{r,c}), continue; end
				[n,m] = size(F{r,c});
				if n == m, fprintf(gp,'\nset size square\n'); end
				fprintf(gp,'set title "%s"\n',ptitle{r,c});
				fprintf(gp,'set xr[-0.5:%d-0.5]\n',m);
				fprintf(gp,'set yr[-0.5:%d-0.5]\n',n);
				fprintf(gp,'set xtics (%s)\n',gp_numtics(1:m,0:m-1));
				fprintf(gp,'set ytics (%s)\n',gp_numtics(1:n,0:n-1));
				fprintf(gp,'set cbrange [0:%g]\n',Fmax(r,c));
				fprintf(gp,'plot "-" matrix with image not\n');
				for i = 1:n
					fprintf(gp,' %12.8f',F{r,c}(i,:)); fprintf(gp,'\n');
				end
				fprintf(gp,'e\ne\n');
			end
		end
		fprintf(gp,'\nunset multiplot\n');
	else
		[n,m] = size(F);
		if n == m, fprintf(gp,'\nset size square\n'); end
		fprintf(gp,'set title "%s"\n',ptitle);
		fprintf(gp,'set xr[-0.5:%d-0.5]\n',m);
		fprintf(gp,'set yr[-0.5:%d-0.5]\n',n);
		fprintf(gp,'set xtics (%s)\n',gp_numtics(1:m,0:m-1));
		fprintf(gp,'set ytics (%s)\n',gp_numtics(1:n,0:n-1));
		fprintf(gp,'set cbrange [0:%g]\n',Fmax);
		fprintf(gp,'plot "-" matrix with image not\n');
		for i = 1:n
			fprintf(gp,' %12.8f',F(i,:)); fprintf(gp,'\n');
		end
		fprintf(gp,'e\ne\n');
	end

	gp_close(gp,gpstem,gpterm);

end

function plot_gc_mat(F,ptitle,cm,Fmax,nr,nc,i)

	if isempty(F), return; end

	if nargin > 4, subplot(nr,nc,i); end

	F = flipud(F);

	[n,m] = size(F);

	colormap(cm);
	imagesc(F,[0 Fmax]);
	axis('square');
	xlabel('from');
	ylabel('to');
	set(gca,'XTick',1:n);
	set(gca,'XTickLabel',1:n);
	set(gca,'YTick',1:m);
	set(gca,'YTickLabel',m:-1:1);
	title(ptitle)

end
