function [g,gsiz,idx] = check_group(group,n)

assert(iscell(group),                     'groups must be specified as a cell vector of indices');
idx = cell2mat(group(:)');
assert(length(unique(idx)) == length(idx),'groups indices must be unique and non-overlapping');
if nargin > 1 && ~isempty(n);
	assert(all(idx >=1 & idx <= n),       'some group indices out of range');
end
g = length(group);
if nargout > 1
	gsiz = cellfun(@length,group);
end
