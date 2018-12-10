function [sparkle_preY] = runSparkleSequence(I,NB,pY,boxSparkle,zoomSparkle,map,beta,WINDOW,sparkleSZ,UNET,UU)
    tm = clock;
    for p = 1:size(pY,1)
        curPoint = (double(round(pY(p,:))));
        BOX = point2Box(curPoint,boxSparkle);
        
        [subI] = mimcrop(I,BOX,zoomSparkle,[]);
        
        
        
        %%%%
        %ssubI = rgb2ind(imresize(subI,.35),map);
        %POS = UNET.predict(ssubI) + UU;
        %%%%
        
        
       
        
        
        POS = fitbitp(subI,NB,map,beta,WINDOW,true);
      
        sparkle_preY(:,:,p) = reshape(POS,sparkleSZ);
        %{
        imshow(subI,[]);
        drawnow
        hold on;
        plot(sparkle_preY(:,1,p),sparkle_preY(:,2,p),'k.')
        hold off
        drawnow
        waitforbuttonpress
        %}
        
        sparkle_preY(:,:,p) = bsxfun(@minus,sparkle_preY(:,:,p),[size(subI,2) size(subI,1)]/2);
        sparkle_preY(:,:,p) = sparkle_preY(:,:,p)*zoomSparkle^-1;
        sparkle_preY(:,:,p) = bsxfun(@plus,sparkle_preY(:,:,p),curPoint);
    end
    etm = etime(clock,tm);
    fprintf(['Done running sparkle sequence in:' num2str(etm) '\n']);
end