function [curSeg] = igetCurveSegment(curve,radiusBundle,radiusSegment)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculates the tangent and normal bundles for the curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:   
    %           curve   : = curve for analysis
    %           radius  : = radius for pca
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           BV      : = a sequence of basis vectors for the tangent and
    %           normal space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp = 1;
    
    numP = size(curve,1) - radiusSegment;
    segP = 2*radiusSegment;
    curSeg = zeros(numP,segP,2);
    
    for i = 1:numP
        tBV = igetFrame_atP(curve,i,radiusBundle);
        
        curPoint = curve(i,:);
        affineT = tBV';
        
        tempCurve = bsxfun(@minus,curve,curPoint);
        
        tempCurve = affineT*tempCurve';
        tempCurve = tempCurve';
    
        dist = sum(tempCurve.*tempCurve,2).^.5;
        
        fidx = find((dist < radiusSegment) & (tempCurve(:,1) >= 0));
        
        
        tmpSeg = tempCurve(fidx,:);
        curSeg(i,:,:) = arcLength(tmpSeg,'spec',segP);
        if disp
            plot(tempCurve(:,1),tempCurve(:,2),'r')
            hold on
            plot(squeeze(curSeg(i,:,1)),squeeze(curSeg(i,:,2)),'b')
            hold off
            axis([-200 200 -200 200])
            drawnow
        end
        
    end
end