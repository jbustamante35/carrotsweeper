function [crv, mline] = getContourAndMidline(msk, vis)
%% getContourAndMidline: extract midline and contour from binary mask image
% This function runs Nathan's algorithm for extracting a contour and calculating
% the midline from a binary mask image.
%
% Usage:
%   [crv, mline] = getContourAndMidline(msk, vis)
%
% Input:
%   msk: binary mask (black object, white background)
%   vis: boolean to visualize output
%
% Output:
%   mline: extracted midline data
%   crv: extracted countour data
%

%% Set constants for respective algorithms
MASK_THRESH = 300; % Original 300
% MIN_THRESH_SIZE = 300; % Remove midlines of certain number of coordinates

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
    %% Some initial processing of the mask
    if size(msk, 3) > 1
        msk = rgb2gray(msk);
    end
    
    msk(:,end) = [];
    msk        = handleFLIP(msk, []);
    msk        = padarray(msk, [0 MASK_THRESH], 'pre', 'replicate');
    skel       = ~msk;
    
    % Run through main functions
    crv   = getBWContour(skel);
    tIdx  = getTipIdx(crv);
    mline = generateMidline(msk, crv, tIdx);
    
    %     % Remove data below threshold number of coordinates [DEPRECATED]
    %     rmIdx          = mline(:,1) < MIN_THRESH_SIZE;
    %     mline(rmIdx,:) = [];
    %     crv(rmIdx,:)   = [];
    %     mline(:,1)     = mline(:,1) - MIN_THRESH_SIZE;
    
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
        imshow(Io, []);
        hold on;
        
        plt(crv.data', 'r.', 5);
        plt(midlines.data', 'g--', 2);
        plt(crv.data(:, tipIdx)', 'g*', 8);
        
    else
        % Subtract midline by size for some reason
        sz         = size(msk);
        mline(:,1) = mline(:,1) - sz(2)/2;
        mline(:,1) = -mline(:,1);
        mline(:,1) = mline(:,1) + sz(2)/2;
        
        % Subtract contour by size for some reason
        crv(:,1) = crv(:,1) - sz(2)/2;
        crv(:,1) = -crv(:,1);
        crv(:,1) = crv(:,1) + sz(2)/2;
        
    end
    
catch e
    fprintf(2, 'Error extracting Midline and Contour\n%s\n', e.getReport);
    mline = [];
    crv  = [];
end

end
