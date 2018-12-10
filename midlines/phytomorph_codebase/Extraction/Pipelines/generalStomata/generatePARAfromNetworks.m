function [para,cP,yP0] = generatePARAfromNetworks(patch,trainedNetwork_Angle_local,Z0)


    yP0 = trainedNetwork_Angle_local.predict(patch);
    yP0 = yP0.*Z0(1,:)+Z0(2,:);
    cP = reshape(yP0(1,1:8),[2 4])';
    para(1:3) = yP0(end-2:end);
    para(4:5) = [0 0];
    %{
    Rpatch = imrotate(patch,-yP0(end)*180/pi,'crop');
    
    yP1 = trainedNetwork_uM.predict(Rpatch);
    yP1 = yP1.*Z1(1,:)+Z1(2,:);

   

   
    %}
end