%% mvgc_makemex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEPRECATED: use Makefile in 'C' subfolder %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build MVGC |mex| files
%
% <matlab:open('mvgc_makemex.m') code>
%
%% Syntax
%
%     mvgc_makemex(force_recompile,verbose)
%
%% Arguments
%
% _input_
%
%     tocompile        name of mex routine to compile (default: 'all')
%     force_recompile  forced recompilation flag (default: false)
%     verbose          verbosity flag (default: false)
%
%% Description
%
% Builds and then tests all MVGC |mex| files from |C| source in the |C|
% subdirectory. If a |mex| file for the current platform already exists in the
% |mex| subdirectory it is just tested, unless the |force_recompile| flag is
% set, in which case it is recompiled. A build is not considered successful
% unless it has also been successfully tested.
%
% Assumes a working Matlab-compatible |C| for your platform which supports |C99|.
% Sensible defaults are hard-coded, but you may have to (or want to) tweak this
% function for your platform/compiler.
%
% _*Note 1:*_ The toolbox is currently distributed with pre-built and tested
% |mex| files for 64-bit Unix (including Linux), Windows and Mac, as these were
% the only test platforms available to us. If Matlab crashes on you, there is
% a very good chance that a pre-built |mex| is to blame. In this case (assuming
% you have a Matlab-compatible C compiler available) you should try running
% <mvgc_makemex.html |mvgc_makemex|> with the |force_recompile| flag set.
%
% _*Note 2:*_ The pre-built Windows 64-bit |mex| files distributed with the
% toolbox were compiled with MinGW64. Note that Microsoft Visual Studio does
% not support the |C99| standard; we reccommend MinGW64,
%
%% See also
%
% <startup.html |startup|> |
% <matlab:doc('mex') |mex|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function mvgc_makemex(tocompile,force_recompile,verbose)

% Build MVGC 'mex' files

if nargin < 1 || isempty(tocompile),       tocompile       = 'all'; end;
if nargin < 2 || isempty(force_recompile), force_recompile = false; end;
if nargin < 3 || isempty(verbose),         verbose         = false; end;

compall = strcmpi(tocompile,'all');

% default mex command

MFLAGS = '-O -largeArrayDims';
if verbose, MFLAGS = ['-v ' MFLAGS]; end % spews lots of details
mexcmd = ['mex ' MFLAGS];

if isunix
    plat = 'Unix';
	CFLAGS = '-std=c99 -march=native -O3 -flto -Wall -Werror -Wextra -Wconversion -Winline -pedantic-errors -D_POSIX_C_SOURCE=199309L -D_BSD_SOURCE -static';
    mexcmd = [mexcmd ' CFLAGS="\$CFLAGS ' CFLAGS '"'];
elseif ispc
    plat = 'Windows';
    CFLAGS = '-std=c99 -Wall -Werror -O3';  % MinGW64
    mexcmd = [mexcmd ' CFLAGS="' CFLAGS '"'];
elseif ismac
    plat = 'Mac';
    % If you want to override compiler flags, you're on your own...
else
    plat = 'Unknown';
    fprintf(2,'\nNOTE: At present ''mvgc_makemex'' has only been tested on Unix (gcc) and Windows (MinGW64).\nIf you are on a different  platform/compiler and have any useful compilation tips, it would be\nvery helpfull if you could tell the MVGC maintainers about it.\n');
end

cc = mex.getCompilerConfigurations();

global mvgc_root;
cdir = fullfile(mvgc_root,'C');
mexdir = fullfile(mvgc_root,'mex');

fprintf('\nYour platform   appears to be : %s (mex extension: %s)\n',plat,mexext);
fprintf('Your C compiler appears to be : %s\n\n',cc(2).Name);

if compall || strcmpi(tocompile,'genvar')

	cfroot   = 'genvar';
	mfroot   = [cfroot '_mex'];
	cfile    = [mfroot '.c'];
	mexfile  = [mfroot '.' mexext];
	mexinvoc = [mexcmd ' ' cdir filesep cfile ' -outdir ' mexdir];

	global have_genvar_mex;
	if exist(mfroot,'file') == 3 && ~force_recompile
		fprintf('\nA mex file ''%s'' already exists for your platform. If you want\nto recompile it, then re-run this routine with the ''force_recompile'' flag set.\n\n',mexfile);
		have_genvar_mex = testmex(cfroot);
	else
		fprintf('\n*** Going to compile ''%s''\n',cfile);
		fprintf('Using mex command ''%s''\n',mexinvoc);
		try
			eval(mexinvoc);
			fprintf('\nLooks like ''%s'' compiled ok. ',cfile);
			have_genvar_mex = testmex(cfroot);
		catch err
			fprintf(2,err.message);
			fprintf(2,'\nHrmmph. ''%s'' failed to compile. Please tell the maintainer about this.\nMeanwhile don''t panic, an equivalent (but slower) scripted ''%s'' will be used.\n\n',cfile,cfroot);
			return
		end
	end

end

if compall || strcmpi(tocompile,'genvma')

	cfroot   = 'genvma';
	mfroot   = [cfroot '_mex'];
	cfile    = [mfroot '.c'];
	mexfile  = [mfroot '.' mexext];
	mexinvoc = [mexcmd ' ' cdir filesep cfile ' -outdir ' mexdir];

	global have_genvma_mex;
	if exist(mfroot,'file') == 3 && ~force_recompile
		fprintf('\nA mex file ''%s'' already exists for your platform. If you want\nto recompile it, then re-run this routine with the ''force_recompile'' flag set.\n\n',mexfile);
		have_genvma_mex = testmex(cfroot);
	else
		fprintf('\n*** Going to compile ''%s''\n',cfile);
		fprintf('Using mex command ''%s''\n',mexinvoc);
		try
			eval(mexinvoc);
			fprintf('\nLooks like ''%s'' compiled ok. ',cfile);
			have_genvma_mex = testmex(cfroot);
		catch err
			fprintf(2,err.message);
			fprintf(2,'\nHrmmph. ''%s'' failed to compile. Please tell the maintainer about this.\nMeanwhile don''t panic, an equivalent (but slower) scripted ''%s'' will be used.\n\n',cfile,cfroot);
			return
		end
	end

end

if compall || strcmpi(tocompile,'findin')

	cfroot   = 'findin';
	mfroot   = [cfroot '_mex'];
	cfile    = [mfroot '.c'];
	mexfile  = [mfroot '.' mexext];
	mexinvoc = [mexcmd ' ' cdir filesep cfile ' -outdir ' mexdir];

	global have_findin_mex;
	if exist(mfroot,'file') == 3 && ~force_recompile
		fprintf('\nA mex file ''%s'' already exists for your platform. If you want\nto recompile it, then re-run this routine with the ''force_recompile'' flag set.\n\n',mexfile);
		have_findin_mex = testmex(cfroot);
	else
		fprintf('\n*** Going to compile ''%s''\n',cfile);
		fprintf('Using mex command ''%s''\n',mexinvoc);
		try
			eval(mexinvoc);
			fprintf('\nLooks like ''%s'' compiled ok. ',cfile);
			have_findin_mex = testmex(cfroot);
		catch err
			fprintf(2,err.message);
			fprintf(2,'\nHrmmph. ''%s'' failed to compile. Please tell the maintainer about this.\nMeanwhile don''t panic, an equivalent (but slower) scripted ''%s'' will be used.\n\n',cfile,cfroot);
			return
		end
	end

end

% More compilations and test functions go here

end

function success = genvar_mex_test

	global have_genvar_mex;
	n = 11;
	p = 23;
	m = 10000;
	rho = 0.95;
	A = specnorm(randn(n,n,p),rho);
	a = squeeze(specnorm(randn(1,1,p),rho));
	oldstate = rng_seed(67132);
	XA = randn(n,m);
	Xa = randn(n,m);
	Xu = randn(1,m);
	rng_restore(oldstate);
	have_genvar_mex = true;
	tic
	YA1 = genvar(A,XA);
	Ya1 = genvar(a,Xa);
	Yu1 = genvar(a,Xu);
	toc
	have_genvar_mex = false;
	tic
	YA2 = genvar(A,XA);
	Ya2 = genvar(a,Xa);
	Yu2 = genvar(a,Xu);
	toc
	maxad = max([maxabs(YA1-YA2) maxabs(Ya1-Ya2) maxabs(Yu1-Yu2)]);
	fprintf('\nMaximum absolute difference = %.g\n',maxad);
	success = maxad < 1e-12;

end

function success = genvma_mex_test

	global have_genvma_mex;
	n = 11;
	p = 23;
	m = 10000;
	rho = 0.95;
	B = specnorm(randn(n,n,p),rho);
	b = squeeze(specnorm(randn(1,1,p),rho));
	oldstate = rng_seed(19676);
	XB = randn(n,m);
	Xb = randn(n,m);
	Xu = randn(1,m);
	rng_restore(oldstate);
	have_genvma_mex = true;
	tic
	YB1 = genvma(B,XB);
	Yb1 = genvma(b,Xb);
	Yu1 = genvma(b,Xu);
	toc
	have_genvma_mex = false;
	tic
	YB2 = genvma(B,XB);
	Yb2 = genvma(b,Xb);
	Yu2 = genvma(b,Xu);
	toc
	maxad = max([maxabs(YB1-YB2) maxabs(Yb1-Yb2) maxabs(Yu1-Yu2)]);
	fprintf('\nMaximum absolute difference = %.g\n',maxad);
	success = maxad < 1e-12;

end

function success = findin_mex_test

	global have_findin_mex;
	a = [8 3 7 1 9];
	b = [6 9 1 3 8 5 7];
	k = [5 4 7 3 2];
	have_findin_mex = true;
	k1 = findin(a,b);
	have_findin_mex = false;
	k2 = findin(a,b);
	assert(isequal(k2,k),'Weird, the scripted version seems to be broken...');
	success = isequal(k1,k);

end

% Test the mex file

function success = testmex(cfroot)

    mfroot  = [cfroot '_mex'];
    mexfile = [mfroot '.' mexext];
    testfun = [cfroot '_mex_test'];

    fprintf('You should test ''%s'' (note: If Matlab crashes the test failed ;-)\n\n',mexfile);
    reply = input('Go for test? y/n [y]: ', 's');

    success = false;
    if isempty(reply) || strcmpi(reply,'y')
        if eval(testfun)
            fprintf('Congratulations, ''%s'' passed the test! It will now be used by default.\n\n',mexfile);
            success = true;
        else
            fprintf(2,'Drat. ''%s'' failed the test. Please tell the maintainer about this.\nMeanwhile don''t panic, an equivalent (but slower) scripted ''%s'' will be used.\n\n',mexfile,cfroot);
        end
    else
        fprintf(2,'\nOk, you bottled out. An equivalent (but slower) scripted ''%s'' will be used.\n\n',cfroot);
    end

end
