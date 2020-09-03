function [sK , pK] = processCurvatures(k, fnm, splt, uids, regx, seps, rgns, sect)
%% processCurvatures: process curvatures into table format
% Description
%
% Usage:
%   [sK , pK] = processCurvatures(k, fnm, regx, seps, rgns, sect)
%
% Input:
%   k: structure of curvatures from all regions and secions
%   fnm: filename from input
%   splt: data is split into upper and lower section
%   uids: metadata to add table columns (default: {UID|Genotype|Replicate})
%   regx: regular expression to separate curvature (default: \s*)
%   seps: character to separate curvatures (default: ';')
%   rgns: regions of curvatures (default: [whole|shoulder|tip])
%   sect: sections of curvatures (default: [upper|lower])
%
% Output:
%   sK: curvatures converted to strings with concatenated fields
%   pK: curvatures converted to strings with fields intact
%

%% Default values
if nargin < 4
    uids = {'UID' , 'Genotype' , 'Replicate'};
    regx = '\s*';
    seps = ';';
    rgns = {'shoulder' , 'tip' , 'whole'};
    sect = {'upper' , 'lower'};
end

%% Convert curvatures to single-line string
if splt
    kstr = @(k,r,s) regexprep(num2str(k.(r).(s)'), regx, seps);
    pK   = cellfun(@(r) cellfun(@(s) kstr(k,r,s), sect, 'UniformOutput', 0)', ...
        rgns, 'UniformOutput', 0)';
    pK   = cell2struct(pK, rgns);
    pK   = structfun(@(r) cell2struct(r, sect), pK, 'UniformOutput', 0);
    
    % Convert to single structure
    S  = cellfun(@(r) cellfun(@(s) kstr(k, r, s), ...
        sect, 'UniformOutput', 0), rgns, 'UniformOutput', 0);
    F  = cellfun(@(r) cellfun(@(s) sprintf('%s_%s', r, s), ...
        sect, 'UniformOutput', 0), rgns, 'UniformOutput', 0);
    
    S = cat(2, S{:});
    F = cat(2, F{:});
else
    kstr = @(k,r,s) regexprep(num2str(k.(r)'), regx, seps);
    pK   = cellfun(@(r) kstr(k,r), rgns, 'UniformOutput', 0)';
    pK   = cell2struct(pK, rgns);
    
    % Convert to single structure
    S  = cellfun(@(r) kstr(k, r), rgns, 'UniformOutput', 0);
    F  = cellfun(@(r) sprintf('%s', r), rgns, 'UniformOutput', 0);
end

%% Extract information from filename and add columns
vals = cellfun(@(x) getNameID({fnm},x), uids, 'UniformOutput', 0);
S    = [vals , S]';
F    = [uids , F]';
sK   = cell2struct(S,F);

end

