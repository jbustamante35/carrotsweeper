function [mline, cntr, sMsk] = cropTraceStraightenSample(msk)
try
    [mline, cntr] = getMidlineAndContour(msk);
    sMsk          = sampleStraighten(mline, flip(msk, 2));
    
catch ME
    getReport(ME)
end
end