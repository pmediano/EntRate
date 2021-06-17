function morder = moselect(stitle,mosel,varargin)

% Report and select model order

ncrits = length(varargin)/2;
assert(iscell(varargin) && isvector(varargin) && 2*(ncrits/2) == ncrits,'Model orders must be spupplied as [''name'',values] pairs');

nlen = zeros(ncrits,1);
for c = 1:2:2*ncrits
    moname = varargin{c};
    moval = varargin{c+1};
    assert(ischar(moname) && isint(moval),'Entry %d: criterion name must be a string and value an integer',c);
    nlen(c) = length(moname);
end
maxnlen = max(nlen);
usersupp = isscalar(mosel) && isint(mosel);
if ~ischar(mosel), assert(usersupp,'Model order selection criterion must be a string or a positive integer'); end
if usersupp, maxnlen = max(maxnlen,4); end % for 'user'

prestr = sprintf('\n\t%%-%ds',maxnlen); % format string for name
morder = [];
fprintf('\n%s:',stitle);
for c = 1:2:2*ncrits
    moname = varargin{c};
    moval = varargin{c+1};
    fprintf([prestr ' = %d'],moname,moval);
    if strcmpi(mosel,moname)
		morder = moval;
		fprintf(' *** selected');
	end
end
if isempty(morder)
    assert(usersupp,'Model order selection criterion ''%s'' doesn''t match any supplied criterion',mosel);
    fprintf([prestr ' = %d *** selected'],'user',mosel);
    morder = mosel;
end
fprintf('\n');
