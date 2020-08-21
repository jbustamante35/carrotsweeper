function [rIdxs , rx] = pullRandom(X, n, getrx)
%% pullRandom: pull random number(s) from distribution
% Description
%
% Usage:
%    rIdxs = pullRandom(X, n, getrx)
%
% Input:
%    X: distribution of numbers
%    n: number of random pulls (optional) [defaults to 1]
%    getrx: return the actual value instead of the index
%
% Output:
%    rIdxs: random index or indices from distribution
%    x: random data pulled from distribution
%
% Author Julian Bustamante <jbustamante@wisc.edu>
%

%% Default to take 1 sample
switch nargin
    case 1
        n     = 1;
        getrx = 0;
    case 2
        getrx = 0;
    case 3
    otherwise
        fprintf(2, 'Error with inputs\n');
        [rIdxs , rx] = deal([]);
        return;
end

%% Return random index and object
% If Shuffle function not available
% rIdxs = randi([1 , length(X)], 1);
% rx    = X(rIdxs);
rIdxs = sort(Shuffle(length(X), 'index', n));
rx    = X(rIdxs);

if getrx
    rIdxs = rx;
end

end


