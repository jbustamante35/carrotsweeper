function [X, map, mask] = readico(filename, index)
%READICO Read icon from a Windows ICO file
%   [X,MAP] = READICO(FILENAME) reads image data from an ICO file
%   containing one or more Microsoft Windows icon resources.  X is
%   a 2-D uint8 array.  MAP is an M-by-3 MATLAB colormap.  If
%   FILENAME contains more than one icon resource, the first will
%   be read.
%
%   [X,MAP] = READICO(FILENAME,INDEX) reads the icon in position INDEX
%   from FILENAME, which contains multiple icons.
%
%   [X,MAP,MASK] = READICO(FILENAME,...) returns the transperency mask
%   for the given image from FILENAME.
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $  $Date: 2011/05/17 02:28:02 $

if (nargin < 2)
    index = 1;
end

if (~ischar(filename))
    error(message('MATLAB:imagesci:readico:badFilename'));
end;

if (isempty(strfind(filename,'.')))
    filename=[filename,'.ico'];
end;

info = imfinfo(filename, 'ico');

if (ischar(index))
    error(message('MATLAB:imagesci:readico:badIndexType', index));
end

if (index < 1)
    error(message('MATLAB:imagesci:readico:badIndex'))
end	

if (index > length(info))
    error(message('MATLAB:imagesci:readico:indexOutOfRange', index, filename, length( info )));
end

% Read the XOR data and its colormap
X = readbmpdata(info(index));

map = info(index).Colormap;

maskinfo = info(index);
maskinfo.BitDepth = 1;

% Calculate the offset of the AND mask.
% Bitmap scanlines are aligned on 4 byte boundaries
imsize = maskinfo.Height * (32 * ceil(maskinfo.Width / 32))/8;
maskinfo.ImageDataOffset = maskinfo.ResourceDataOffset + ...
    maskinfo.ResourceSize - imsize;

mask = readbmpdata(maskinfo);
