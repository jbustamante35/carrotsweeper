function [rIdxs , val] = pullRandom(X, n, getval)
%% pullRandom: pull random number(s) from distribution
% Description
%
% Usage:
%    [rIdxs , val] = pullRandom(X, n, getval)
%
% Input:
%    X: distribution of numbers
%    n: number of random pulls (optional) [defaults to 1]
%    getval: return the actual value instead of the index
%
% Output:
%    rIdxs: random index or indices from distribution
%    val: random data pulled from distribution
%
% Author Julian Bustamante <jbustamante@wisc.edu>
%

%% Default to take 1 sample
switch nargin
    case 1
        n      = 1;
        getval = 0;
    case 2
        getval = 0;
    case 3
    otherwise
        fprintf(2, 'Error with inputs\n');
        [rIdxs , val] = deal([]);
        return;
end

%% Return random index and object
% Check if Shuffle function is available
if ~isempty(which('Shuffle'))
    rIdxs = sort(Shuffle(length(X), 'index', n));
else
    % No Shuffle function found
    rIdxs = sort(randi(length(X), [1 , n]));
end
val = X(rIdxs);

if getval
    rIdxs = val;
end

end
