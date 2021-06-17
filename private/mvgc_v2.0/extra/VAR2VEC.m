function [PI,G] = VAR2VEC(A)

PI = -(eye(size(A,1))-sum(A,3));

if nargout > 1
    p = size(A,3);
    G = -A(:,:,2:p);
    for k = p-2:-1:1
        G(:,:,k) = G(:,:,k)+G(:,:,k+1);
    end
end
