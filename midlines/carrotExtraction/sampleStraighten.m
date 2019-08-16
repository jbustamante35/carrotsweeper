function strVec = sampleStraighten(mline, msk)
%% sampleStraighten: midline-based straightener
% This function straightens an object in an image by extending from the midline
% coordinates and extracting the pixel values each coordinate corresponds to.
%
% Usage:
%   strVec = sampleStraighten(mline, msk, crv)
%
% Input:
%   mline: midline coordinates
%   msk: binary mask image
%   crv: contour to define boundary to omit pixels during straightening
%
% Output:
%   strVec: straightened vector structure
%

try
    %%
    fmsk     = flip(msk, 3);
    [dS, dG] = extendCarrotMidline(mline, [0 0], fmsk);
    dSize    = size(dG);
    chnls    = size(msk, 3);
    
    if chnls > 1
        strVec     = zeros(length(dS), chnls);
        
        %% Interpolation of mask pixels along midline
        % Multi-channel image
        for chnl = 1 : chnls
            strVec(:,chnl) = ...
                ba_interp2(double(msk(:,:,chnl)) / 255, dS(:,2), dS(:,1));
        end

    else
        % Binary image
        mmsk   = double(fmsk) / 255;
        strVec = ba_interp2(mmsk, dS(:,2), dS(:,1));
    end
    
    %% Reshape to generate straight mask
    strVec = imbinarize(reshape(strVec, dSize(1:2))', ...
        'adaptive', 'Sensitivity', 1);    
    
    %% Run from midline, extend from normal and remove out-of-contour pixels

    
catch e
    fprintf(2, 'Error straightening mask\n%s\n', e.getReport);
    strVec = [];
    
end

end