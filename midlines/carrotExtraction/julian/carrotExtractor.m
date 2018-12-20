function [mline, cntr, smsk, pmsk] = carrotExtractor(dataIn, vis, svData, svFigs)
%% carrotExtractor: midline extraction and straightener
% This is a detailed description of this script...
%
% Note about saving with the svData parameter:
% When svData is set to true, the output is saved in a .mat file that is placed
% in the output directory [ see Usage below ]. All 4 output variables are stored 
% in a single structure called CARROTS. 
%
% Usage:
%   [mline, cntr, smsk, pmsk] = carrotExtractor(dataIn, vis, svData, svFigs)
%
% Input:
%   dataIn: path to directory of binary images
%   vis: boolean visualize outputs
%   svData: boolean to save output in .mat file [ see note above ]
%   svFigs: boolean to save figures of straightened carrots
%
% Output:
%   mline: cell array of midline data
%   cntr; cell array of contour data
%   smsk; cell array of straightened mask images
%   pmsk: cell array of processed mask images
%
% Usage (continued):
% If you want to save the results, the data will automatically be placed in the
% input directory as a subfolder named output_yymmdd, where 'yymmdd' corresponds 
% to the year (y), month (m), and today's date (d).
%
% dataIn = '/path/to/directory/of/masks';
%

%% Load file list of binary mask images
% dataIn = '~/Dropbox/EdgarSpalding/projects/carrotsweeper/data/development_181011';

if svData || svFigs
    dOut    = sprintf('output_%s', tdate('s'));
    dataOut = sprintf('%s/%s', dataIn, dOut);
    mkdir(dataOut);
end

ext = '.png';
fin = imageDatastore(dataIn, 'FileExtensions', ext);

%% Extract Midline, Contour, Straightened Image, Straightened Mask
tot                       = numel(fin.Files);
[mline, cntr, pmsk, smsk] = deal(cell(1, tot));
for n = 1 : tot
    try
        % Prepare mask for extraction functions
        pmsk{n} = extendDimension(fin.readimage(n), 0);
        pmsk{n} = double(imcomplement(pmsk{n}));
        
        % Run processed mask through extraction functions
        [mline{n}, cntr{n}] = getMidlineAndContour(pmsk{n}, vis);
        smsk{n}             = sampleStraighten(mline{n}, flip(pmsk{n}, 2), vis);
        
        % Clear figure axis
        if vis && n < tot
            cla;clf;
        end
        
    catch e
        fprintf(2, 'Error in Carrot Extraction Pipeline\n%s\n', e.getReport);
    end
end

%% Show Output of processed and straightened masks
if vis
    psec  = 0.7;
    fig   = figure;
    set(fig, 'Color', 'w');
    
    for n = 1 : tot
        try            
            plotCarrots(n, pmsk{n}, mline{n}, cntr{n}, smsk{n}, psec);
            
        catch e
            fprintf(2, 'Error plotting figure for data %d\n%s\n', ...
                n, e.getReport);
            
            plotCarrots(n, pmsk{n}, [0 0], [0 0], [0 0], psec);
            
        end
        
        % Save figures in output directory
        if svFigs
            fName = getDirName(dataIn);
            fnm   = sprintf('%s/straightCarrot%d_%s', dataOut, n, fName);
            savefig(fig, fnm);
            saveas(fig, fnm, 'tiffn');
        end
        
        cla;clf;
    end
end

%% Save Data in output directory
if svData
    fName   = getDirName(dataIn);
    CARROTS = v2struct(mline, cntr, smsk, pmsk);
    nm      = sprintf('%s/%s_carrotExtractor_%s_%dCarrots', ...
        dataOut, tdate('s'), fName, tot);
    save(nm, '-v7.3', 'CARROTS');
end

end

function pad = extendDimension(msk, val)
%% Extend dimension of mask by defined number of pixels in val
% The threshold parameter THERSH defines theminimum length of either dimension.
% If a dimension is lower than the threshold value, it is padded to match the
% threshold value with the value defined in the val parameter. This function
% automatically detects which dimension should be padded.
THRESH = 200;
chk    = size(msk) - THRESH;
idx    = chk < 0;

if idx
    pad = msk;
else
    pad = padarray(msk, abs(chk(idx)), val);
end

end

function plotCarrots(idx, raw_mask, midline, contours, straight_mask, psec)
%% plotCarrots: plotting function for this script
subplot(211);
img = handleFLIP(raw_mask, 3);
imshow(img, []);
hold on;
plt(midline, 'r-', 2);
plt(contours, 'b-', 2);
ttlP = sprintf('Midline and Contour on Mask\nCarrot %d', idx);
title(ttlP);

subplot(212);
flp = handleFLIP(straight_mask, 2);
imshow(flp, []);
ttlS = sprintf('Straighted Mask\nCarrot %d', idx);
title(ttlS);

pause(psec);

end