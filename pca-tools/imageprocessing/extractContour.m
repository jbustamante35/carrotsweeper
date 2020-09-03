function [cout , cntr] = extractContour(bw, ncrds, nrm, alt, req)
%% extractContour: find contour of single image
% This function blah
%
% Usage:
%   [cout , cntr] = extractContour(bw, ncrds, nrm, alt, req)
%
% Input:
%   bw: bw image
%   ncrds: number of coordinates to normalize contour
%   nrm: reindexing method for NormalizedOutline
%   alt: use alternative reindexing parameter to use HypoQuantyl's method
%   req: string to output '', 'Interp', or 'Normalized' Outline
%
% Output:
%   cout: just the outline, not the ContourJB object
%   cntr: various data from contour at given frame
%

%% Handle inputs
switch nargin
    case 1
        ncrds = 800;
        nrm   = 'default';
        alt   = 'default';
        req   = 'Interp';
    case 2
        nrm   = 'default';
        alt   = 'default';
        req   = 'Interp';
    case 3
        alt   = 'default';
        req   = 'Interp';
    case 4
        req   = 'Interp';
    case 5
    otherwise
        fprintf(2, 'Error with inputs [Expected 5 | Given %d]\n', nargin);
        [cntr, cout] = deal([]);
        return;
end

%% Get boundaries of inputted bw image
bndAll   = bwboundaries(bw, 'noholes');
[~, lrg] = max(cellfun(@numel, bndAll));
bnds     = bndAll{lrg};

%% Interpolate distances to an equalized number of coordinates
bnds   = [getDim(bnds, 2) , getDim(bnds, 1)]; % Switch y-/x-coordinates to x-/y-coordinates
intrps = interpolateOutline(bnds, ncrds);

%% Output final structure
% Normalization method after reindexing: 'default' subtracts by mean after
% reindexing, while 'alt' just reindexes
switch alt
    case 'default'
        % Default algorithm [CarrotSweeper] with requested normalization method
        cntr = ContourJB('Outline', bnds, 'InterpOutline', intrps);
        cntr.ReindexCoordinates(nrm);
    case 'alt'
        % Alternative start coordinate by adding alt parameter [HypoQuantyl]
        cntr = ContourJB('Outline', bnds, 'InterpOutline', intrps, 'AltInit', alt);
        cntr.ReindexCoordinates(nrm);
    otherwise
        fprintf(2, ...
            'Error re-indexing coordinates. Must be [default|alt] (%s)\n', alt);
        [cntr , cout] = deal([]);
        return;
end

%% Ouput the contour requested individually
cname = sprintf('%sOutline', req);
cout  = cntr.(cname);
end


