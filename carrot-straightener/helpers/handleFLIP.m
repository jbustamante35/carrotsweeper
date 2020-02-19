function [img , direc2flip] = handleFLIP(img, direc2flip)
%% handleFLIP
% Description
%
% Usage:
%    [img , direc2flip] = handleFLIP(img, direc2flip)
%
% Input:
%    img: image to flip
%    direc2flip: direction (1,2,3,4) to flip [empty to auto-detect]
%
% Output:
%    img: flipped image
%    direc2flip: direction image was flipped
%
% Author Nathan Miller <nbmill@gmail.com>
% Edited by Julian Bustamante <jbustamante@wisc.edu>
%

%%
if isempty(direc2flip)
    % handle the direction
    [direc2flip] = DC(img);
end

switch direc2flip
   case 1                                                      % from the left
       % do nothing
   case 2                                                      % from the top
       img = img';
   case 3                                                      % from the right
       img = fliplr(img);
   case 4                                                      % from the bottom
       img = flipud(img)'; % WORKS
end
% handle the direction
% read and flip
