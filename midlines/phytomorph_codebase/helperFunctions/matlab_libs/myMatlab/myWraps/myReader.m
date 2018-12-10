function [I] = myReader(varargin)
    if nargin >= 1;file=varargin{1};end
    if nargin >= 2;para=varargin{2};end
    if nargin == 1
        % read image
        I = imread(file);
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % set variables
        cenP = para{1};
        half_width = para{2};
        half_height = para{3};
        COL = round([(cenP(1) - half_width) (cenP(1) + half_width)]);
        ROW = round([(cenP(2) - half_height) (cenP(2) + half_height)]);
        % set variables
        I = double(imread(file,'PixelRegion',{ROW,COL}));
    end
end