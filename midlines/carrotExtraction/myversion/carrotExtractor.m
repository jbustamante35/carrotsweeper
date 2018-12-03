function [mline, cntr, smsk, pmsk] = carrotExtractor(vis, svData, svFigs)
%% carrotExtractor: midline extraction and straightener
%
%
% Input:
%   vis: boolean visualize outputs
%   svData: boolean to save output in .mat file
%   svFigs: boolean to save figures of straightened carrots
%
% Output:
%   mline: cell array of midline data
%   cntr; cell array of contour data
%   smsk; cell array of straightened mask images
%   pmsk: cell array of processed mask images
%

%% Load file list of binary mask images
fldr   = '/home/jbustamante/Dropbox/EdgarSpalding/projects/carrotsweeper';
din    = 'data/development_181011/binary_mask';
dataIn = sprintf('%s/%s', fldr, din);

if svData || svFigs
    dOut    = sprintf('output_%s', tdate('s'));
    dataOut = sprintf('%s/%s', dataIn, dOut);
    mkdir(dataOut);
end

ext  = '.png';
fin  = imageDatastore(dataIn, 'FileExtensions', ext);

%% Extract Midline, Contour, Straightened Image, Straightened Mask
tot                       = numel(fin.Files);
[mline, cntr, pmsk, smsk] = deal(cell(1, tot));
for n = 1 : tot
    try
        % Prepare mask for extraction functions
        pmsk{n}             = extendDimension(fin.readimage(n), 0);
        pmsk{n}             = double(imcomplement(pmsk{n}));
        
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
        subplot(211);
        imshow(pmsk{n});
        ttlP = sprintf('Processed Mask\nCarrot %d', n);
        title(ttlP);
        
        subplot(212);
        flp = handleFLIP(smsk{n}, 2);
        imshow(flp);
        ttlS = sprintf('Straighted Mask\nCarrot %d', n);
        title(ttlS);
        
        pause(psec);
        
        % Save figures in output directory
        if svFigs
            fName = getDirName(din);
            fnm   = sprintf('%s/straightCarrot%d_%s', dataOut, n, fName);
            savefig(fig, fnm);
            saveas(fig, fnm, 'tiffn');
        end
    end
end

%% Save Data in output directory
if svData
    fName   = getDirName(din);
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