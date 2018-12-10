function [isUrl, filenameOut] = getFileFromURL(filenameIn)
%GETFILEFROMURL Detects whether the input filename is a URL and downloads
%file from the URL

%   Copyright 2007-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/05/17 02:27:20 $

% Download remote file.
if (strfind(filenameIn, '://'))
  
    isUrl = true;

    if (~usejava('jvm'))
        error(message('MATLAB:imagesci:getFileFromURL:noJVM'))
    end
    
    try
        filenameOut = urlwrite(filenameIn, tempname);
    catch
        error(message('MATLAB:imagesci:getFileFromURL:urlRead', filenameIn));
    end
    
else
  
    isUrl = false;
    filenameOut = filenameIn;
    
end
