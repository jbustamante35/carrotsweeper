function [crv , cntr] = getBWContour(msk, CNTR_LENGTH)
%% getBWContour: extracts a ContourJB from a binary mask image
% Wrapper of my method to extract contours from binary mask images. Output is in
% the form of Nathan's curve structure and my ContourJB object.
%
% Usage:
%   [crv , cntr] = getBWContour(msk, CNTR_LENGTH)
%
% Input:
%   msk: binary masks image
%   CNTR_LENGTH: interpolation size for number of coordinates per contour
%
% Output:
%   crv: x-/y-coordinates of the extracted contour
%   cntr: my custom ContourJB class
%

% Interpolation size for extracted contour
if nargin < 2; CNTR_LENGTH = 800; end

[crv , cntr] = extractContour(msk, CNTR_LENGTH);
end