function [ X ] = SimulateEvenProcess(T)
%%SIMULATEEVENPROCESS Generate data from R. James' even process
%
%   X = SIMULATEEVENPROCESS(T) simulates a time series of length T from the
%   even process as defined by James et al, 2011. The even process has an
%   entropy rate of 2/3 bit.
%
%
% Reference:
%   James, Ellison & Crutchfield (2011). Anatomy of a bit: Information in a
%   time series observation. Chaos 21, 037109.
%
% Pedro Mediano and Fernando Rosas, Jun 2021

X = zeros([1, T]);
X(1) = 1*(rand() < 0.5);
counter = 0;
for i=2:T

%   if X(i-1) && mod(counter, 2) == 1
%     X(i) = true;
%     counter = counter + 1;
%   else
%     X(i) = rand() < 0.5;
%   end

    if X(i-1) == 1
        counter = counter + 1;
    elseif X(i-1) == 0
        counter = 0;
    end

    if mod(counter, 2) == 0
        X(i) = 1*(rand() < 0.5);
    else
        X(i) = 1;
    end

end

X = logical(X);

