function [endPoints5 , endPoints10] = testSkeletonize(mline, msk)
% testSkeletonize: 
% Scott's implementation to convert coordinates to logical, then perform
% skeletonization and de-spurring on midline coordinates.
%
% Input:
%    mline: midline coordinates
%    msk: binary mask associated with midline
%
% Output:
%    endPoints5: Coordinates of midline's end points
%    endPoints10: What's the difference from endPoints5?
%

%% Save midlines 
% mline10 doesn't branch, mline5 does
figure
fig = plt(mline{5}, 'k-', 5);
axis off;
saveas(fig,'mline5','png');

figure
fig = plt(mline{10}, 'k-', 5);
axis off;
saveas(fig,'mline10','png');

%% Read them in, convert to grayscale, invert 1s and 0s, and binarize
image5    = imread('mline5.png');
image5    = rgb2gray(image5);
image5    = imcomplement(image5);
BW5       = imbinarize(image5);

image10   = imread('mline10.png');
image10   = rgb2gray(image10);
image10   = imcomplement(image10);
BW10      = imbinarize(image10);

%% Skeletonize and de-spur little sideshoots
spurCount = 20;
BW5       = bwmorph(BW5,'skel',inf);
BW5       = bwmorph(BW5, 'spur',spurCount);
BW10      = bwmorph(BW10,'skel',inf);
BW10      = bwmorph(BW10, 'spur',spurCount);

%% Count the number of endpoints
endPoints5 = bwmorph(BW5, 'endpoints');
sum(endPoints5(:) == 1) %%3

endPoints10 = bwmorph(BW10, 'endpoints');
sum(endPoints10(:) == 1) %%2

end