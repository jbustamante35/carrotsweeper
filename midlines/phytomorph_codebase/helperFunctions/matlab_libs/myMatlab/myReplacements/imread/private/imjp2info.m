function metadata = imjp2info(filename)
%IMJP2NFO Information about a JPEG 200 file.
%   METADATA = IMJP2INFO(FILENAME) returns a structure containing
%   information about the JPEG 2000 file specified by the string
%   FILENAME.  
%
%   See also IMFINFO.

%   Copyright 2008-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/17 02:27:48 $

% JPEG2000 is not supported on Solaris.	 
if (isequal(computer(), 'SOL64'))	 
    error(message('MATLAB:imagesci:imjp2info:unsupportedPlatform'))	 
end	 

% Call the interface to the Kakadu library.
metadata = imjp2infoc(filename);

d = dir(filename);
metadata.Filename = filename;
metadata.FileModDate = datestr(d.datenum);
metadata.FileSize = d.bytes;

