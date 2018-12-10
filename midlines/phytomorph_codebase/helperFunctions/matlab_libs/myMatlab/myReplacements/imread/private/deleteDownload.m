function deleteDownload(filename)
%DELETEDOWNLOAD Deletes the temporary file downloaded from the URL

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/05/17 02:27:19 $

try
    delete(filename);
catch
    warning(message('MATLAB:imagesci:deleteDownload:removeTempFile', filename))
end
