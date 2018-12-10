function [phDocument] = grapeClusterMetric_ver1(Mask)
    Mask = bwlarge(Mask);
    R = regionprops(Mask,'Area');
    massDistributionVer = sum(Mask,2);
    massDistributionHor = sum(Mask,1);
    phDocument.Area = R(1).Area;
    phDocument.massDistributionVer = massDistributionVer;
    phDocument.massDistributionHor = massDistributionHor;
end
