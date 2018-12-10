function [updatedPara,displaceMent,newPara] = suggestNewPara(currentPara,cP,patchCenter,w)
    
   

    majorAxis = cP(1,:)-cP(3,:);
    minorAxis = cP(2,:)-cP(4,:);

    centerPointMajor = .5*majorAxis;
    centerPointMinor = .5*minorAxis;
    displaceMent = mean(cat(1,centerPointMinor,centerPointMajor),1);

    gM = norm(majorAxis)/2;
    gu = norm(minorAxis)/2;
    
    esti(1,:) = cP(1,:) - displaceMent;
    esti(2,:) = cP(2,:) - displaceMent;
    esti(3,:) = -(cP(3,:) - displaceMent);
    esti(4,:) = -(cP(4,:) - displaceMent);


    ang(1) = atan2(esti(1,2),esti(1,1));
    ang(2) = atan2(esti(2,2),esti(2,1))-pi/2;
    ang(3) = atan2(esti(3,2),esti(3,1));
    ang(4) = atan2(esti(4,2),esti(4,1))-pi/2;

    for l = 1:size(esti,1)
        quiver(displaceMent(1),displaceMent(2),esti(l,1),esti(l,2))
hold on
    end
   


    newPara = [mean(ang) gM gu currentPara(4:5)];
    [newPara] = spinOutAngle(newPara);
    %newPara(1) = newPara(1) + currentPara(1);
    updatedPara = newPara*(w) + currentPara*(1-w);
    %updatedPara = mean(cat(1,newPara,currentPara),1);
    

end