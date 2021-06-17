function bedit(fstr)

editor = getenv('MATLAB_EDITOR');

s = which(fstr);

if isempty(s)
    fprintf(2,'Sorry, didn''t find ''%s''\n',fstr);
    return
end

[status res] = system([editor ' ' s ' &']);
if status ~= 0
    fprintf(2,'ERROR(%d): %s\n',status,res);
end
