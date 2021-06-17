function A = var_rand(n,p,rho,w,plotm)

% Generate random VAR covariance matrices
%
% n       - observation variable dimension, or connectivity matrix/array
% p       - number of lags
% rho     - spectral radius
% w       - decay weighting parameter: empty (default) = don't weight
%
% A       - VAR coefficients array
%
% TODO

% plotm = []      - don't plot
% plotm = n       - Matlab plot to figure n (if zero, use next)
% plotm = string  - Gnuplot terminal (may be empty)

if nargin < 4, w     = []; end
if nargin < 5, plotm = []; end

if isscalar(n)
    A = randn(n,n,p);
else
    C = n; % connectivity matrix
    [n,n1,q] = size(C);
    assert(n1 == n,'Connectivity matrix must be square');
    if q == 1
		C = repmat(C,[1 1 p]);
	else
		assert(isempty(p),'If full connectivity array, model order must be empty');
		p = q;
	end
    A = C.*randn(n,n,p);
end

if isempty(w)
	A = specnorm(A,rho);
else
	A = specnorm(exp(-w*sqrt(p))*A,rho);
end

if ~isempty(plotm) % we're going to plot
    a = zeros(p,1);
    for k = 1:p
        a(k) = norm(A(:,:,k));
    end
    if ischar(plotm) % plotm is gpterm
        gp_qplot((1:p)',a,[],'set xlabel "k"\nset ylabel "|A|"\nset xtics 1\nunset key\nset grid',plotm);
    elseif plotm
		if plotm == 0, figure; else, figure(plotm); end; clf;
        plot((1:p)',a);
        ylim([0 1]);
        xlabel('$k$','Interpreter','latex');
        y=ylabel('$\Vert A_k \Vert$','rot',0,'Interpreter','latex');
        set(y, 'position', get(y,'position')-[0.2,0,0]);
    end
end
