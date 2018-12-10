function [M,CM] = readxwd(fname)
%READXWD Read image data from an XWD file.
%   [X,MAP] = READPCX(FILENAME) reads image data from a PCX file.
%   X is a 2-D uint8 array.  MAP is an M-by-3 MATLAB-style
%   colormap.  XWDREAD can the following types of XWD files:
%    1-bit or 8-bit ZPixmaps
%    XYBitmaps
%    1-bit XYPixmaps
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/05/17 02:28:10 $

%   Reference:  Murray and vanRyper, Encyclopedia of Graphics
%   File Formats, 2nd ed, O'Reilly, 1996.

info = imxwdinfo(fname);

if (info.ByteOrder == 0)
    [fid, fopenMessage] = fopen(fname, 'r', 'ieee-le');
elseif (info.ByteOrder == 1)
    [fid, fopenMessage] = fopen(fname, 'r', 'ieee-be');
else
    error(message('MATLAB:imagesci:readxwd:badByteOrder'));
end
if (fid == -1)
    error(message('MATLAB:imagesci:readxwd:fileOpen', fname, fopenMessage));
end

% See if it is an XYPixmap (unsupported)
if (strcmp(info.PixmapFormat, 'XYPixmap') && (info.PixmapDepth ~= 1))
    error(message('MATLAB:imagesci:readxwd:badXYPixmap', info.PixmapDepth));
end

fseek(fid, info.HeaderSize, 'bof');
tmp = fread(fid,[6 info.NumColormapEntries],'uint16');
CM = [(tmp(2,:)+tmp(1,:)*65535)',tmp(3:5,:)'/65535];

% What storage format is the data in?
if (strcmp(info.PixmapFormat, 'ZPixmap') || ...
            strcmp(info.PixmapFormat, 'XYBitmap') || ...
            (strcmp(info.PixmapFormat, 'XYPixmap') && (info.PixmapDepth == 1)))
    % ZPixmap,  XYbitmap or XYpixmap depth 1
    % Size of pixels
    if (ismember(info.PixmapDepth, [1 8]))
        prec = 'uint8';
    else
        error(message('MATLAB:imagesci:readxwd:badZPixmap', info.PixmapDepth));
    end

    if (info.PixmapDepth == 8), % ZPixmap with depth 8
        pad = ceil(info.Width/info.BitmapPad*info.PixmapDepth) * ...
                info.BitmapPad/info.PixmapDepth - info.Width;
        M = uint8(fread(fid,[info.Width+pad info.Height], prec))';
    else      
        % ZBitmap XYPitmap or XYPixmap of depth 1
        M = fread(fid,[info.Width/8 info.Height],'uint8')';

        % Take apart bytes into 8 bits and store as a matrix of 1's and 0's
        [a,b]=size(M);
        P=logical(repmat(uint8(0), a, b*8));
        for i=7:-1:0,
            tmp=find(M/2^i >= 1);
            if ~isempty(tmp),
                M(tmp)=M(tmp)-2^i;
                P(tmp+ceil(tmp/a)*(7*a)+(-i)*a+a*(2*i-7) * ...
                        (~info.BitmapBitOrder)) = 1;
            end
        end
        CM=[0 1 1 1;1 0 0 0];
        M=P;
    end
    
else
    error(message('MATLAB:imagesci:readxwd:badFormat', info.PixmapDepth, info.PixmapFormat));
end

CM=CM(:,2:4);

fclose(fid);
