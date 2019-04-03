function [crv, cntr] = getBWContour(msk)
%% getBWContour: extracts a ContourJB from a binary mask image
% Wrapper of my method to extract contours from binary mask images. Output is in
% the form of Nathan's curve structure and my ContourJB object.
%
% Usage:
%   [crv, cntr] = getBWContour(msk)
%
% Input:
%   msk: binary masks image
%
% Output:
%   crv: x-/y-coordinates of the extracted contour
%   cntr: my custom ContourJB class
%

%%
CNTR_LENGTH = 800; % Interpolation size for extracted contour
cntr        = extractContour(msk, CNTR_LENGTH);
crv         = cntr.InterpOutline;

end