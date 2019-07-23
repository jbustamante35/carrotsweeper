function [skel, crv, mline, tCrds, dsts] = getContourAndMidline(msk, vis)
%% getContourAndMidline: extract midline and contour from binary mask image
% This function runs Nathan's algorithm for extracting a contour and calculating
% the midline from a binary mask image.
%
% NOTE: msk input should be black object on white background!
% Usage:
%   [skel, crv, mline] = getContourAndMidline(msk, vis)
%
% Input:
%   msk: binary mask (must be black object, white background)
%   vis: boolean to visualize output
%
% Output:
%   skel: binary mask flipped to face left-right
%   crv: extracted countour data
%   mline: extracted midline data
%   tCrds: coordinates of the tip
%   dsts: vector of distance transform along midline
%

%% Set constants for respective algorithms
MASK_THRESH     = 100; % Original 300
MIN_THRESH_SIZE = 100; % Remove midlines of certain number of coordinates
FACE            = 3;   % Direction to point images (original 3)

% TODO: Pass these constants to the respective algorithm
% getBWContour
%     CNTR_LENGTH = 800; % Interpolation size for extracted contour
%
% getTipIdx
%     KSNIP   = 163;
%     SMOOTH1 = 5;
%     SMOOTH2 = 10;
%
% generateMidline
%     DSK   = 31;
%     FSPEC = 7;
%

try
    %% Initial processing of the mask to face left-right and then pad
    if size(msk, 3) > 1
        msk = rgb2gray(msk);
    end
    
    % Force flip to left-right [arg = 3]
    %     skel = handleFLIP(msk, []);
    skel = handleFLIP(msk, FACE);
    chk  = ~imbinarize(skel);
    
    % Remove rows that are all empty
%     cols2remove = 1;
    cols2remove = 5;
    while ~sum(chk(:,1))
        chk(:,1)    = [];
        cols2remove = cols2remove  + 1;
    end
    
    skel = skel(:, cols2remove:end);
    skel = ~padarray(skel, [0 MASK_THRESH], 'pre', 'replicate');
    
    %% Run through main functions
    crv = getBWContour(skel);
    
    % Remove padded area of mask and curve/midline coordinates
    rmCrv        = crv(:,1) < MIN_THRESH_SIZE;
    crv(rmCrv,:) = [];
    
    % Identify tip as point of highest curvature
    %     tCrds = getTipIdx(skel);
    tCrds = getTipIdx(chk);
    tCrds = [tCrds(:,1) + MIN_THRESH_SIZE tCrds(:,2)];
    
    % Generate midline starting from tip and distance transform values
    [mline, dsts] = generateMidline(~skel, tCrds);
    
    %% Remove padding from mask, contour, and midline
    % Remove padded area of mask
    skel(:,1:MIN_THRESH_SIZE) = [];
    crv(:,1)                  = crv(:,1) - MIN_THRESH_SIZE;
    
    % Remove any midline coordinates beyond mask
    rmMid          = mline(:,1) < MIN_THRESH_SIZE;
    mline(rmMid,:) = [];
    mline(:,1)     = mline(:,1) - MIN_THRESH_SIZE;
    
    % Shift tip coordinates to new location
    tCrds = [tCrds(:,1) - MIN_THRESH_SIZE , tCrds(:,2)];
    
    %% Visualize output
    if vis
        cla;clf;
        imagesc(skel);
        colormap gray;
        axis image;
        hold on;
        
        plt(crv, 'b.', 4);
        plt(mline, 'r.', 3);
        plt(tCrds, 'g*', 8);
        %         bar(flip(dsts), 1, 'r');
    end
    
catch e
    fprintf(2, 'Error extracting Midline and Contour\n%s\n', e.getReport);
    mline = [];
    crv   = [];
end

end
