function [cropLine] = getCarrotSplitLine(kidx,blueClusterNumber)
    splitLine = kidx == blueClusterNumber;
    splitLine = bwlarge(splitLine);
    sig = sum(splitLine,1);
    cropLine = find(sig > 10);
    cropLine = cropLine(1);
end