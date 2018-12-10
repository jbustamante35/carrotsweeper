function list = parselist(listin)
%Parse comma separated list into a cell array

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2010/08/07 07:31:49 $

if (isempty(listin))
    list = {};
else
    % Return a row vector.
    list = textscan(listin, '%s', 'delimiter', ',');
    list = list{1}';
end
