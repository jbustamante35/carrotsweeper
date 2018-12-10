function [imgOut] = subtractImages2( foreground, background);
% SUBTRACTIMAGES2 reads the foreground and background images and returns
% the foreground subtracted from the background.
    %   Detailed Explanation

fore = double(imread(foreground))/255;
back = double(imread(background))/255;

imgOut = back - fore;
end