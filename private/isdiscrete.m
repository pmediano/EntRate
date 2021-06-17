function [ D ] = isdiscrete(X)
%%ISDISCRETE True for arrays of discrete variable-compatible types.
%
%   ISDISCRETE(X) returns true if X is an array of a data type suitable to be
%   interpreted by ELPh as a discrete random variable. This is true for integer,
%   logical, and categorical types, as well as for float and double type arrays
%   in which all numbers are close to integers.
%
% Examples:
%   isdiscrete([true false]);  % returns true
%
%   isdiscrete(int8(2));       % returns true
%
%   isdiscrete([3 4 5]);       % returns true
%                              % though note that isinteger([3 4 5]) returns
%                              false, because numbers are stored as double by
%                              default
%
% Pedro Mediano, Feb 2021

if ischar(X) || isstring(X) || isstruct(X)
  error('ELPh:DataType', ['String, char, and struct data not accepted. ', ...
    'Please transform your data to double or int.']);
elseif ~exist('OCTAVE_VERSION', 'builtin') && iscategorical(X)
  D = true;
else
  D = isinteger(X) || islogical(X) || sum(abs(X(:) - round(X(:)))) < 1e-10;
end

