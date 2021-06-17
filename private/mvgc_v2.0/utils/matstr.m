function s = matstr(A,fmtstr)

% Return a numeric matrix as a formatted string

assert(ismatrix(A));
s = sprintf([repmat(fmtstr,1,size(A,2)) '\n'],A);
