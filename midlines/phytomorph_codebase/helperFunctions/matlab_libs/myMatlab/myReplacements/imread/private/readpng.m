function [X,map,alpha] = readpng(filename, varargin)
%READPNG Read an image from a PNG file.
%   [X,MAP] = READPNG(FILENAME) reads the image from the
%   specified file.
%
%   [X,MAP] = READPNG(FILENAME,'BackgroundColor',BG) uses the
%   specified background color for compositing transparent
%   pixels.  By default, READPNG uses the background color
%   specified in the file, if present.  If not present, the
%   default is either the first colormap color or black.  If the
%   file contains an indexed image, BG must be an integer in the
%   range [1,P] where P is the colormap length.  If the file
%   contains a grayscale image, BG must be an integer in the
%   range [0,65535].  If the file contains an RGB image, BG must
%   be a 3-element vector whose values are in the range
%   [0,65535].
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2011/05/17 02:28:06 $

bg = parse_args(varargin{:});

if (isempty(bg) && (nargout >= 3))
    % User asked for alpha and didn't specify a background
    % color; in this case we don't perform the compositing.
    bg = 'none';
end

alpha = [];
[X,map,h] = pngreadc(filename, bg);
X = permute(X, ndims(X):-1:1);

% A one-row image needs to have this dimension re-imposed.
if h == 1
	X = reshape(X,[1 size(X)]);
end

if (ismember(size(X,3), [2 4]))
    alpha = X(:,:,end);
    % Strip the alpha channel off of X.
    X = X(:,:,1:end-1);
end



%--------------------------------------------------------------------------
function bg = parse_args(param,value)
bg = [];

if nargin < 1
    return
end

% Process param/value pairs.  Only 'backgroundcolor' is recognized.
if (~ischar(param))
    error(message('MATLAB:imagesci:readpng:parameterNotString'));
end

n = numel(param);
if ~strncmpi(param,'backgroundcolor',n)
    error(message('MATLAB:imagesci:readpng:unrecognizedParameter', param));
end

bg = value;
return

