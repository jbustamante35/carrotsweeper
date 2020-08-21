function vals = getNameID(FNAMES, id)
%% getNameID: extract values from an id in a filename
%
% Usage:
%   vals = getNameID(FNAMES, id)
%
% Input:
%   FNAMES: filename string array
%   id: id string to search for in filename
%
% Output:
%   vals: values of id string from filename
%

expr = sprintf('%s_(?<id>.*?)}', id);
val  = regexpi(FNAMES, expr, 'names');
vals = cellfun(@(x) char(x.id), val, 'UniformOutput', 0);

end