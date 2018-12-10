function [] = cornerExtractor(I)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find all corner(s) in the first image
    cornersFeatures = cornerFeatureMap(5);
    cornerFeatures = cornersFeatures.computeFeatureMap(double(I));
    allCorners = simplePointExtractor(20,.000008);
    allCorners = allCorners.extractPoints(cornerFeatures);
    
    
    
    for e = 1:numel(allCorners)
        pt = allCorners{e};
        
    end
end