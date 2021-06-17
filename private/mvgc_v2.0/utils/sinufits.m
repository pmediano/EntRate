% SINUSOIDAL LEAST-SQUARES-FIT SIGNAL
%
% INVOCATION
%
% [s,p,q] = sinufits(x,fs,f)
%
% INPUTS
%
% x         matrix of time-series values (observations x variables)
% f         sinusoidal frequency to remove (in range (0 .. Nyqvist])
% fs        sampling frequency (default: angular frequency in range 0 ... 2*pi)
%
% OUTPUTS
%
% s         the sinusoidal least-squares fit signal
% p,q       sinusoidal coefficients: fitted sinusoid is p.*cos(w*t) + q.*sin(w*t)
%           NOTE: sinusoid amplitude is hypot(p,q)
%
% Author: Lionel Barnett, University of Sussex, UK.
%
% Date:   November, 2016
%
% lionelb@sussex.ac.uk
% http://users.sussex.ac.uk/~lionelb/
%
% You may do as you choose with this code - a credit would be nice, though :)
%
%%

function [s,p,q] = sinufits(x,fs,f)

if isempty(fs), fs = 2*pi; end % If no sampling frequency supplied assume angular frequency.

assert(ismatrix(x) && isnumeric(x),'Time series must be a numeric matrix');

[T,n] = size(x);

assert(f > eps && f <= fs/2,'Frequency must lie in range (0 .. Nyqvist]');

w  = 2*pi*f/fs; % convert to angular frequency
mu = mean(x);   % save mean
x  = x-mu;      % temporal de-mean
t  = (0:T-1)';  % time sequence

% Calculate the sinusoidal fit

W    = (sin(w*T)/sin(w))/T;
u    = W*cos(w*(T-1));
v    = W*sin(w*(T-1));
comt = cos(w*t);
somt = sin(w*t);
xc   = mean(x.*comt);
xs   = mean(x.*somt);
D    = 2/(1-W*W);
p    = D*((1-u)*xc-v*xs); % cos coefficient
q    = D*((1+u)*xs-v*xc); % sin coefficient
s    = p.*comt+q.*somt;   % the sinusoidal signal
