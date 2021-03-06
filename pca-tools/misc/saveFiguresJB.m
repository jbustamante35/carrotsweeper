function saveFiguresJB(figs, fnms, sav_fig, img_type, sav_dir)
%% saveFiguresJB: save figures as figure names in .fig and .tiffn format
% It's about time I made this function instead of copy-pasting the same simple
% for loop every day.
%
% Usage:
%   saveFiguresJB(figs, fnms, sav_fig, img_type, sav_dir)
%
% Input:
%   figs: numeric array of figure handles
%   fnms: cell array of figure names (should equal number of figs)
%   sav_fig: save .fig format (default to true)
%   img_type: image format to save figures (default to 'tiffn')
%   sav_dir: directory to save figures and images
%

%% Save them
% Default to tiffn format (uncompressed tiff)
switch nargin
    case 1
        fnms     = arrayfun(@(x) sprintf('%s_figure%d', tdate, x), ...
            figs, 'UniformOutput', 0);
        sav_fig  = 0;
        img_type = 'tiffn';
        sav_dir  = pwd;
    case 2
        sav_fig  = 0;
        img_type = 'tiffn';
        sav_dir  = pwd;
    case 3
        img_type = 'tiffn';
        sav_dir  = pwd;
    case 4
        sav_dir  = pwd;
    case 5
    otherwise
        fprintf(2, 'Too many input arguments [%d]\n', nargin);
        return;
end

%% Save figures (create directory if it doesn't exist)
if ~isfolder(sav_dir)
    mkdir(sav_dir);
end

fnms = cellfun(@(fnm) sprintf('%s%s%s', sav_dir, filesep, fnm), ...
    fnms, 'UniformOutput', 0);

if sav_fig
    arrayfun(@(fig) savefig(figs(fig), fnms{fig}), ...
        1 : numel(figs), 'UniformOutput', 0);
    arrayfun(@(fig) saveas(figs(fig), fnms{fig}, img_type), ...
        1 : numel(figs), 'UniformOutput', 0);
else
    arrayfun(@(fig) saveas(figs(fig), fnms{fig}, img_type), ...
        1 : numel(figs), 'UniformOutput', 0);
end

end

