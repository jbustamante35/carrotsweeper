function [] = displayPhenotypes(oI,Mask,phDocument,oPath,baseName)
    
    image(double(oI)/255);
    axis off
    hold on
    plot(phDocument.massDistributionVer,1:numel(phDocument.massDistributionVer),'r');
    plot(1:numel(phDocument.massDistributionHor),phDocument.massDistributionHor,'r');
    hold off
    
    fileName = [oPath baseName '_massDistribution.jpg'];
    saveas(gca,fileName);
    pause(.2);
    
    
    close all
    
    out = flattenMaskOverlay(double(oI)/255,Mask);
    image(out)
    axis off
    hold on
    plot(phDocument.massDistributionVer,1:numel(phDocument.massDistributionVer),'r');
    plot(1:numel(phDocument.massDistributionHor),phDocument.massDistributionHor,'r');
    hold off
    
    fileName = [oPath baseName '_massDistribution_MaskOverlay.jpg'];
    saveas(gca,fileName);
    pause(.2);
    
end