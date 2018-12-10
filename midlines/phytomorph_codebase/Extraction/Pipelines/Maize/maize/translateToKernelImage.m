function [fileName] = translateToKernelImage(plateName,wellName)
    p = '/mnt/spaldingdata/kernelImages/Pictures/';
    plateName = strrep(plateName,'-','_');
    if strcmp(plateName(end),'A') | strcmp(plateName(end),'B')
        plateName(end-1) = [];
    end
    fileName = [p plateName filesep 'bw' filesep plateName '_' wellName(2) wellName(4) '_T.tiff'];
end