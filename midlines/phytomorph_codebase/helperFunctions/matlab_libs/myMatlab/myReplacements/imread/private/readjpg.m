function [A,junk] = readjpg(filename)
%READJPG Read image data from a JPEG file.
%   A = READJPG(FILENAME) reads image data from a JPEG file.
%   A is a uint8 array that is 2-D for grayscale and 2-D for RGB
%   images. 
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Steven L. Eddins, June 1996
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/05/17 02:28:04 $

depth = jpeg_depth(filename);

if (depth <= 8)
  
    A = rjpg8c(filename);
    
elseif (depth <= 12)
  
    A = rjpg12c(filename);
    
elseif (depth <= 16)
  
    A = rjpg16c(filename);
    
else
  
    error(message('MATLAB:imagesci:readjpg:unsupportedJPEGBitDepth', depth))
    
end

junk = [];
