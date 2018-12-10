function info = imrasinfo(filename)
%IMRASINFO Get information about the image in a RAS file.
%
%   INFO = IMRASINFO(FILENAME) returns information about the image
%   contained in a RAS file.
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   A complete official specification for the RAS (Sun Raster) image format
%   does not seem to have been made publicly available.	 As sources for the
%   RAS image format I have used
%
%     * /usr/include/rasterfile.h of the Sun OS
%     * The rasterfile(4) man page of the Sun OS
%     * The files libpnm4.c, rasttopnm.c, and pnmtorast.c in the NetPBM 8.3
%	distribution
%     * "Inside SUN Rasterfile", a note by Jamie Zawinski
%      <jwz@teak.berkeley.edu> containing an excerpt from "Sun-Spots
%      Digest", Volume 6, Issue 84.
%
% Author:	  Peter J. Acklam
% E-mail:	  pjacklam@online.no

%  Copyright 2001-2010 The MathWorks, Inc.
%  $Revision: 1.1.6.6 $  $Date: 2011/05/17 02:27:53 $

% Try to open the file for reading.
fid = fopen(filename, 'rb', 'ieee-be');
if (fid < 0)
    error(message('MATLAB:imagesci:imrasinfo:fileOpen', filename));
end

% Initialize universal structure fields to fix the order
info = initializeMetadataStruct('RAS', fid);

info.FormatSignature = [];

% Initialize RAS-specific structure fields to fix the order.
info.Length	  = [];
info.Type	  = [];
info.MapType	  = [];
info.MapLength = [];

% magic number
info.FormatSignature = fread(fid, 1, 'uint32');

% width (pixels) of image
info.Width		= fread(fid, 1, 'uint32');

% height (pixels) of image
info.Height		= fread(fid, 1, 'uint32');

% depth (1, 8, 24, or 32 bits) pr pixel
info.BitDepth	= fread(fid, 1, 'uint32');

% length (in bytes) of image
info.Length		= fread(fid, 1, 'uint32');

% type of file; see READRAS for details.
info.Type		= fread(fid, 1, 'uint32');

% type of colormap; see READRAS for details.
info.MapType		= fread(fid, 1, 'uint32');

% length (bytes) of following map
info.MapLength	= fread(fid, 1, 'uint32');

% see if we got the whole header
if isempty(info.MapLength)
    fclose(fid);
    error(message('MATLAB:imagesci:imrasinfo:truncatedHeader')); 
end

% we have got what we need from the file, so close it
fclose(fid);

% get the color type
if info.MapLength
    info.ColorType = 'indexed';
elseif info.BitDepth == 1
    info.ColorType = 'grayscale';
else
    info.ColorType = 'truecolor';
end
