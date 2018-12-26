function [mline, cntr] = getMidlineAndContour(msk, vis)
%% getMidlineAndContour: extract midline and contour from binary mask image
% This function runs Nathan's algorithm for extracting a contour and calculating 
% the midline from a binary mask image. 
% 
% Usage:
%   [mline, cntr] = getMidlineAndContour(msk, vis)
%
% Input:
%   msk: binary mask (black object, white background)
%   vis: boolean to visualize output
%
% Output:
%   mline: extracted midline data
%   cntr: extracted countour data
% 

%%
try
    % 
    out = isolate_carrot_Roots(msk, 0, [], []);    
    
    %
    mline       = out(1).midlines.data';    
    rm          = mline(:,1) < 300;
    mline(rm,:) = [];
    mline(:,1)  = mline(:,1) - 300;
    
    %
    cntr       = out(1).contours.data';
    rm         = cntr(:,1) < 300;
    cntr(rm,:) = [];
    cntr(:,1)  = cntr(:,1) - 300;
    
    %% Visualize output
    if ~vis
        % Subtract midline by size for some reason
        sz         = size(msk);
        mline(:,1) = mline(:,1) - sz(2)/2;
        mline(:,1) = -mline(:,1);
        mline(:,1) = mline(:,1) + sz(2)/2;
        
        % Subtract contour by size for some reason
        cntr(:,1) = cntr(:,1) - sz(2)/2;
        cntr(:,1) = -cntr(:,1);
        cntr(:,1) = cntr(:,1) + sz(2)/2;
    end
catch e
    fprintf(2, 'Error extracting Midline and Contour\n%s\n', e.getReport);
    mline = [];
    cntr  = [];
end

end
