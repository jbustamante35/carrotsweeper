function info = impnginfo(filename)
%IMPNGNFO Information about a PNG file.
%   INFO = IMPNGINFO(FILENAME) returns a structure containing
%   information about the PNG file specified by the string
%   FILENAME.  
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Steven L. Eddins, August 1996
%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $  $Date: 2011/05/17 02:27:51 $

try
    info = pnginfoc(filename);
    s = dir(filename);
    info.FileModDate = datestr(s.datenum);
    info.FileSize = s.bytes;
catch myException
    error(message('MATLAB:imagesci:impnginfo:libraryFailure', myException.message));
end

