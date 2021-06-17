function mextest

% Test MVGC 'mex' files

if isunix
    plat = 'Unix';
elseif ispc
    plat = 'Windows';
elseif ismac
    plat = 'Mac';
else
    plat = 'Unknown';
end

cc = mex.getCompilerConfigurations();

global mvgc_root;
mexdir = fullfile(mvgc_root,'mex');

global have_mvfilter_mex; have_mvfilter_mex = false;
global have_findin_mex;   have_findin_mex   = false;

mexdir = fullfile(mvgc_root,'mex');

fprintf('\nYour platform   appears to be : %s (mex extension: %s)\n',plat,mexext);
fprintf('Your C compiler appears to be : %s\n\n',cc(2).Name);

testmex('mvfilter',mexdir);
testmex('findin',  mexdir);

end

function testmex(cfroot,mexdir)

	mfroot  = [cfroot '_mex'];
	mexfile = [mfroot '.' mexext];
	testfun = [cfroot '_mex_test'];

	fprintf('*** Testing ''%s'' (note: If Matlab crashes the test failed ;-)\n\n',mexfile);

	if exist(fullfile(mexdir,mexfile),'file') ~= 3
		fprintf(2,'ERROR: mex file ''%s'' not found.\n\n',mexfile);
		return
	end

	if eval(testfun)
		fprintf('Congratulations: ''%s'' passed the test! It will now be used by default.\n\n',mexfile);
	else
		fprintf(2,'ERROR: ''%s'' failed the test. Please tell the maintainer about this.\nMeanwhile don''t panic, an equivalent (but slower) scripted ''%s'' will be used.\n\n',mexfile,cfroot);
	end
end

function success = mvfilter_mex_test

	global have_mvfilter_mex;
	n = 11;
	p = 23;
	q = 15;
	m = 10000;
	rhoa = 0.95;
	rhob = 0.85;
	A = specnorm(randn(n,n,p),rhoa);
	a = squeeze(specnorm(randn(1,1,p),rhoa));
	B = specnorm(randn(n,n,q),rhob);
	b = squeeze(specnorm(randn(1,1,q),rhob));
	oldstate = rng_seed(67132);
	X = randn(n,m);
	rng_restore(oldstate);
	have_mvfilter_mex = true;
	fprintf('mex version      : '); tic
	YX1 = mvfilter(B,A,X);
	YX2 = mvfilter(b,A,X);
	YX3 = mvfilter(B,a,X);
	YX4 = mvfilter(b,a,X);
	YX5 = mvfilter([],A,X);
	YX6 = mvfilter([],a,X);
	YX7 = mvfilter(B,[],X);
	YX8 = mvfilter(b,[],X);
	t = toc; fprintf('%8.6f seconds\n',t);
	have_mvfilter_mex = false;
	fprintf('scripted version : '); tic
	YM1 = mvfilter(B,A,X);
	YM2 = mvfilter(b,A,X);
	YM3 = mvfilter(B,a,X);
	YM4 = mvfilter(b,a,X);
	YM5 = mvfilter([],A,X);
	YM6 = mvfilter([],a,X);
	YM7 = mvfilter(B,[],X);
	YM8 = mvfilter(b,[],X);
	t = toc; fprintf('%8.6f seconds\n',t);
	maxad = max([maxabs(YX1-YM1) maxabs(YX2-YM2) maxabs(YX3-YM3) maxabs(YX4-YM4) maxabs(YX5-YM5) maxabs(YX6-YM6) maxabs(YX7-YM7) maxabs(YX8-YM8)]);
	fprintf('\nMaximum absolute difference = %.g\n',maxad);
	success = maxad < 1e-12;
	if success, have_mvfilter_mex = true; end

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
	if success, have_findin_mex = true; end

end
