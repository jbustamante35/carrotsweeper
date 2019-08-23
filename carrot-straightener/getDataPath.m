function [fout, dout] = getDataPath(droot, d)
%% getDataPath: get paths of output data from a root directory
% Date string (parameter d) should be in format outputted by my tdate function.
% So it should be a character array of 'yymmdd' (i.e. 190611)
%
% Input:
%   droot: directory for single genotype
%   d: date string to search through

outdir = sprintf('output-%s', d);
outtyp = '.mat';

% Get output directory
lst      = dir2(droot);
outdirs  = lst(cell2mat(arrayfun(@(x) x.isdir, lst, 'UniformOutput', 0)));
datadir  = outdirs(cell2mat(arrayfun(@(x) contains(x.name, outdir), ...
    outdirs, 'UniformOutput', 0)));
dout     = sprintf('%s/%s', datadir.folder, datadir.name);

% Get .mat file from output directory
fins    = dir(dout);
datafin = fins(cell2mat(arrayfun(@(x) contains(x.name, outtyp), ...
    fins, 'UniformOutput', 0)));
fout    = sprintf('%s/%s', datafin.folder, datafin.name);

end