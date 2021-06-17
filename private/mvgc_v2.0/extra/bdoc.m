function bdoc(fstr)

browser = getenv('MATLAB_BROWSER');

if nargin < 1 || isempty(fstr)
    dpath = fullfile(matlabroot,'help','matlab','index.html'); % display Matlab help index
    system([browser ' ' dpath ' &']);
    return
end

docroot = fullfile(getenv('LOCALREPO'),'matlab','ncomp','mvgc','docs','html');  % first look in MVGC docs
if disdoc(docroot,fstr,browser), return; end

docroot = fullfile(matlabroot,'help','matlab','ref'); % then look in 'matlab/ref' section
if disdoc(docroot,fstr,browser), return; end

docroot = fullfile(matlabroot,'help','matlab'); % then look in 'matlab' section
if disdoc(docroot,fstr,browser), return; end

docroot = fullfile(matlabroot,'help'); % then look in entire 'help' section
if disdoc(docroot,fstr,browser), return; end

fprintf(2,'Sorry, didn''t find any help for ''%s''\n',fstr);

function success = disdoc(docroot,fstr,browser)

[status, dpath] = system(['find ' docroot ' -iname ' fstr '.html']);
assert(status == 0,'''find'' call failed');
success = ~isempty(dpath);
if success
    dpaths = strsplit(dpath);
    dpath = dpaths{1};
    %web(['file://' dpaths{1}],'-browser');
    system([browser ' ' dpath ' &']);
end
