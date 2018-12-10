function out = impixel_figtools(img)
%IMPIXEL_FIGTOOLS Pixel color values with figure tools on.
%   IMPIXEL_FIGTOOLS works very much like the MATLAB's in-built function
%   IMPIXEL and lets you use the important figure toolbar features like
%   Zoom In, Zoom Out, Pan, Rotate and Data Cursor, which are not
%   available with IMPXEL. As of now, IMPIXEL_FIGTOOLS only supports
%   functionality of inputting one input, which has to be the image data.
%
%   In the syntax below, IMPIXEL_FIGTOOLS displays the input image and
%   waits for you to specify the pixels with the mouse:
%   P = IMPIXEL_FIGTOOLS(I)
%   Here P is the matrix of pixel values and I is the image data.
%
%   When you finish selecting pixels, IMPIXEL_FIGTOOLS returns an M-by-3
%   matrix of RGB values in the supplied output argument. If you
%   do not supply an output argument, IMPIXEL_FIGTOOLS returns the matrix
%   in ANS.
%
%   Class Support
%   -------------
%   The input image can be uint8. Though not tested, I would assume it to
%   work well with these formats too - uint16, int16, double, single and
%   logical.
%
%   The output P is double.
%
%   Example
%   -------
%       RGB = imread('peppers.png');
%       pixels = impixel_figtools(RGB)
%
%   See also IMPIXEL.
%
%   Feedback / Bugs / How-this-helped / How-this-sucked /
%   How-this-could_be_improved /Anything about it are MOST welcome.
%
%   Platform: MATLAB R2011B
%
%   Divakar Roy   2012

out = double(cell2mat(ginput_gui(img)));

return;
