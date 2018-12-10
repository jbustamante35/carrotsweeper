function wjpgc(varargin)
%WJPGC Write image data to a JPEG file.
%   WJPGC(IM,FILENAME) writes the grayscale or rgb
%   data in the uint8 array IM to the JPEG file specified
%   by the string FILENAME.
%
%   WJPGC(IM,FILENAME,QUALITY) uses the specified
%   quality factor.
%
%   See also IMWRITE, IMREAD, and IMFINFO.

%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/02/26 17:32:04 $
%#mex

error('MATLAB:imagesci:wjpgc:missingMEX', 'Missing MEX-file WJPGC');
