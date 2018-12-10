function info = imtifinfo(filename)
%IMTIFINFO Information about a TIFF file.
%   INFO = IMTIFINFO(FILENAME) returns a structure containing
%   information about the TIFF file specified by the string
%   FILENAME.  If the TIFF file contains more than one image,
%   INFO will be a structure array; each element of INFO contains
%   information about one image in the TIFF file.
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.12 $  $Date: 2011/05/17 02:27:54 $

% This function is a MATLAB program that does *not* use a libtiff-based
% MEX-file because the im*info functions are used when trying to
% guess the format of an input file because the user didn't
% specify the format.  We don't want to take the memory hit of
% loading a big MEX-file just to see if we have a TIFF file.

% TIFF files might be little-endian or big-endian.  Start with
% little-endian.  If we're wrong, we'll catch it down below and
% reopen the file.
[fid, msg] = fopen(filename, 'r', 'ieee-le');
if (fid == -1)
    error(message('MATLAB:imagesci:imtifinfo:fileOpen', filename, msg));
end



%
% Check that it is a valid tiff file.
sig = fread(fid, 4, 'uint8')';
if (isequal(sig, [73 73 43 0]) || isequal(sig, [77 77 0 43]))
	% We do not as of yet handle Big Tiff
    fclose(fid);
    error(message('MATLAB:imagesci:imtifinfo:bigTiffNotSupported'));
end
if (~isequal(sig, [73 73 42 0]) && ...
    ~isequal(sig, [77 77 0 42]))
    fclose(fid);
    error(message('MATLAB:imagesci:imtifinfo:notTIFF'));
end

if (sig(1) == 73)
    byteOrder = 'little-endian';
else
    byteOrder = 'big-endian';
end


fclose(fid);


raw_tags = tifftagsread ( filename  );
if numel(raw_tags) == 0
    error(message('MATLAB:imagesci:imtifinfo:noImages'));
end
info = tifftagsprocess ( raw_tags );
