% Return indices of first occurrences of elements of row vector a in row
% vector b, in row vector k. Assumes elements are integers, since double
% float comparisons are exact. Returned indices are 1-offset, or zero if
% element of a not found.
%
% NOTE: No error checking at all! Call it right!

function k = findin(a,b)

global have_findin_mex;

if have_findin_mex
	k = findin_mex(a,b); % fast and dirty (for smallish vectors)
else
	[~,k] = ismember(a,b);
end
