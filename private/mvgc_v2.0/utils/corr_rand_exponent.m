function [vexp,iters] = corr_rand_exponent(n,gtarget,iscorr,S,tol,imax,verb)

	if nargin < 3 || isempty(iscorr), iscorr = false; end
	if nargin < 4 || isempty(S),      S      = 1000;  end
	if nargin < 5 || isempty(tol),    tol    = 1e-05; end
	if nargin < 6 || isempty(imax),   imax   = 100;   end
	if nargin < 7 || isempty(verb),   verb   = false; end

	if iscorr
		gtarget = -(n/2)*log(1-gtarget*gtarget);
	end

	g = zeros(S,1);

	pmax = 0.5;
	gmean = 0;
	k = 1;
	while gmean < gtarget
		pmax = 2*pmax;
		gmean = meang(n,pmax,S);
		if verb, fprintf('init %2d : pmax = %7.4f, g = %8.4f\n',k,pmax,gmean); end
		k = k+1;
	end
	if verb, fprintf('\n'); end

	pmin = 0;
	vexp = Inf;
	dvexp = Inf;
	iters = 1;
	while iters <= imax && abs(dvexp) > tol
		ovexp = vexp;
		vexp = (pmin+pmax)/2;
		dvexp = vexp-ovexp;
		gmean = meang(n,vexp,S);
		if verb, fprintf('trial %3d : vexp = %7.4f, g = %8.4f, dvexp = % g\n',iters,vexp,gmean,dvexp); end
		if gmean < gtarget
			pmin = vexp;
		else
			pmax = vexp;
		end
		iters = iters+1;
	end

end

function gmean = meang(n,vexp,S)

	for s = 1:S
		[Q,R] = qr(randn(n));
		M = Q*diag(sign(diag(R))); % M orthogonal
		v = realpow(abs(randn(n,1)),vexp);
		g(s) = sum(log(diag(M*diag(v)*M'))-log(v));
	end
	gmean = mean(g);
end
