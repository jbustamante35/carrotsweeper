function [cropLine] = getCarrotSplitLine(kidx,blueClusterNumber)
    splitLine = kidx == blueClusterNumber;
    splitLine = bwlarge(splitLine,4);
    splitLine = imclose(splitLine,strel('square',101));
    splitLine = bwlarge(splitLine,1);
    sig = sum(splitLine,1);
    cropLine = find(sig > 10);
    cropLine = cropLine(1);
end