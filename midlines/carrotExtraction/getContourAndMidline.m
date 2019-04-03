function [skel, crv, mline] = getContourAndMidline(msk, vis)
%% getContourAndMidline: extract midline and contour from binary mask image
% This function runs Nathan's algorithm for extracting a contour and calculating
% the midline from a binary mask image.
%
% Usage:
%   [skel, crv, mline] = getContourAndMidline(msk, vis)
%
% Input:
%   msk: binary mask (black object, white background)
%   vis: boolean to visualize output
%
% Output:
%   skel: binary mask flipped to face left-right
%   crv: extracted countour data
%   mline: extracted midline data
%

%% Set constants for respective algorithms
MASK_THRESH = 300; % Original 300
MIN_THRESH_SIZE = 300; % Remove midlines of certain number of coordinates

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
    
    skel = handleFLIP(msk, []);
    skel = ~padarray(skel, [0 MASK_THRESH], 'pre', 'replicate');
    
    %% Run through main functions
    crv   = getBWContour(skel);
    tCrds = getTipIdx(crv);
    mline = generateMidline(~skel, tCrds)';
        
    %% Remove padded area of mask and curve/midline coordinates
    rmCrv          = crv(:,1) < MIN_THRESH_SIZE;
    crv(rmCrv,:)   = [];
    rmMid          = mline(:,1) < MIN_THRESH_SIZE;
    mline(rmMid,:) = [];
    
    % Remove padded area of mask
    skel(:,1:MIN_THRESH_SIZE) = [];
    mline(:,1)                = mline(:,1) - MIN_THRESH_SIZE;
    crv(:,1)                  = crv(:,1) - MIN_THRESH_SIZE;
    
    %% Original method [DEPRECATED]
    %     % Run old contour extractor and midline generator
    %     out   = isolate_carrot_Roots(msk, 0, 0);
    %
    %     % Remove output with small output
    %     MIN_THRESH_SIZE = 300; % Remove midlines of certain number of coordinates
    %     crv       = out(1).contours.data';
    %     rm         = crv(:,1) < MIN_THRESH_SIZE;
    %     crv(rm,:) = [];
    %     crv(:,1)  = crv(:,1) - MIN_THRESH_SIZE;
    %
    %     % Remove output with small output
    %     mline       = out(1).midlines.data';
    %     rm          = mline(:,1) < MIN_THRESH_SIZE;
    %     mline(rm,:) = [];
    %     mline(:,1)  = mline(:,1) - MIN_THRESH_SIZE;
    
    %% Visualize output
    if vis
        cla;clf;
        imshow(skel, []);
        hold on;
        
        plt(crv, 'b.', 4);
        plt(mline, 'r.', 3);
        plt(tCrds, 'g*', 8);
        
    else
%         % Subtract midline by size for some reason
%         sz         = size(msk);
%         mline(:,1) = mline(:,1) - sz(2)/2;
%         mline(:,1) = -mline(:,1);
%         mline(:,1) = mline(:,1) + sz(2)/2;
%         
%         % Subtract contour by size for some reason
%         crv(:,1) = crv(:,1) - sz(2)/2;
%         crv(:,1) = -crv(:,1);
%         crv(:,1) = crv(:,1) + sz(2)/2;
        
    end
    
catch e
    fprintf(2, 'Error extracting Midline and Contour\n%s\n', e.getReport);
    mline = [];
    crv  = [];
end

end
