function Mb = memavail

% Take with a pinch of salt

if isunix
	[~,mem] = system('free -m | grep Mem');
	stats = str2double(regexp(mem,'[0-9]*','match'));
	% This is MemAvailable from /proc/meminfo. Maybe MemFree is more approriate?
	% Problem is we don't know the size of largest available contiguous block, or
	% if that even makes sense for Matlab
	Mb = stats(6);
elseif ispc
	user = memory;
	% Okay, but maybe MemAvailableAllArrays is what you want?
	Mb = user.MaxPossibleArrayBytes/1e6;
elseif ismac
	fprintf(2,'SORRY: ''memavail'' not available on MacOS\n');
	Mb = NaN;
else
	error('Unknown operating system');
end
