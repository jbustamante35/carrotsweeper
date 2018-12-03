function [mline, cntr, sMsk, pMsk] = carrotExtractor(vis)
%% Load file list of binary mask images
fldr  = '/home/jbustamante/Dropbox/EdgarSpalding/projects/carrotsweeper';
din   = 'data/development_181011/binary_mask';
dOut    = 'output';

dataIn = sprintf('%s/%s', fldr, din);
dataOut = sprintf('%s/%s', dataIn, dOut);
ext  = '.png';
fin  = imageDatastore(dataIn, 'FileExtensions', ext);

%% Extract Midline, Contour, Straightened Image, Straightened Mask
[mline, cntr, pMsk, sMsk] = deal(cell(1, numel(fin.Files)));
for e = 1 : numel(fin.Files)
    try
        pMsk{e}             = extendDimension(fin.readimage(e), 0);
        pMsk{e}             = double(imcomplement(pMsk{e}));
        [mline{e}, cntr{e}] = getMidlineAndContour(pMsk{e}, vis);
        sMsk{e}             = sampleStraighten(mline{e}, flip(pMsk{e}, 2), vis);
        
        if vis && e < numel(fin.Files)
            cla;clf;
        end
        
    catch ME
        fprintf(2, 'Error\n%s\n', ME.getReport);
    end
end


end

function pad = extendDimension(msk, val)
%% Extend dimension of mask by defined number of pixels
THRESH = 100;
chk = size(msk) - THRESH;
idx = chk < 0;

if idx
    pad = msk;
else
    pad = padarray(msk, chk(idx), val);
end

end