function [X, junk] = readfits(filename)
%READFITS Read image data from a FITS file.
%   A = READFITS(FILENAME) reads the unscaled data from the primary HDU
%   of a FITS file.
%
%   See also FITSREAD.

%   Copyright 1984-2010 The MathWorks, Inc. 
%   $Revision: 1.1.6.4 $  $Date: 2011/05/17 02:27:59 $


warning(message('MATLAB:imagesci:readfits:use_fitsread'));

X = fitsread(filename, 'raw');
junk = [];
