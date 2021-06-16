function p = plt(varargin)
%% plt: my general plotting function for data with x-/y-coordinates
% I've needed this shortcut enough to justify writing a full function for it. User defines the data
% containing [x y] coordinates as main parameter. Optional are the marker type, marker color, and
% marker size. Default values for marker are the default for the plot() function.
%
% Usage:
%   p = plt(D, typ, sz)
%
% Input:
%   D: 2D data containing x-/y-coordinate
%   typ: marker color and/or marker type (optional)
%   sz: marker size (optional)
%
% Output:
%   p: figure handle for plotted data

switch nargin
    case 1
        D = force2D(varargin{1});
        p = plot(D(:,1), D(:,2));
        
    case 2
        D   = force2D(varargin{1});
        typ = varargin{2};
        p   = plot(D(:,1), D(:,2), typ);
        
    case 3
        D = force2D(varargin{1});
        typ = varargin{2};
        sz  = varargin{3};
        
        if contains(typ, '-')
            if numel(sz) > 1
                msz = sz(1);
                lsz = sz(2);
                p = plot(D(:,1), D(:,2), typ, ...
                    'LineWidth', lsz, 'MarkerSize', msz);
            else
                p = plot(D(:,1), D(:,2), typ, 'LineWidth', sz);
            end
        else
            p = plot(D(:,1), D(:,2), typ, 'MarkerSize', sz);
        end
        
    otherwise
        fprintf(2, 'Incorrect number of input\n');
        return;
end

end

function D = force2D(D)
%% force2D: force 1D data to be 2D

if size(D,2) == 1
    D(:,2) = 1:size(D,1);
    D      = fliplr(D);
end

end
