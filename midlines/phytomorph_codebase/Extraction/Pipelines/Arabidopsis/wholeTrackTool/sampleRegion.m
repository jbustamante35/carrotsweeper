function [Z Z1] = sampleRegion(MAP_E,MAP_L,MSK,Q)
Z = zeros(size(MAP_E,1),size(MAP_E,2),size(MAP_E,4));
Z1 = zeros(size(MAP_E,1),size(MAP_E,2),size(MAP_E,4));
[p2,p1] = find(MSK);
idx = find(MSK==1);
for e = 1:size(MAP_E,4)
    LAG = [];
    
    for d = 1:size(MAP_L,3)
        LAG(:,:,d) = imfilter(MAP_L(:,:,d,e),fspecial('gaussian',[31 31],11),'replicate');
    end
    
    
    
    
    np1 = ba_interp2(LAG(:,:,1),p1,p2);
    np2 = ba_interp2(LAG(:,:,2),p1,p2);
    
    
    qn1 = ba_interp2(LAG(:,:,1),Q(1),Q(2));
    qn2 = ba_interp2(LAG(:,:,2),Q(1),Q(2));
    
    
    EUR = [];
    for d = 1:size(MAP_E,3)
        EUR(:,:,d) = imfilter(MAP_E(:,:,d,e),fspecial('gaussian',[31 31],11),'replicate');
    end
    
    
    
    Q1 = ba_interp2(EUR(:,:,1),qn1,qn2);
    Q2 = ba_interp2(EUR(:,:,2),qn1,qn2);
    
  
    
    EUR_SPEED = sum(EUR.^2,3).^.5;
    
    QS = ba_interp2(EUR_SPEED(:,:),qn1,qn2);
    
    EUR_Q_SPEEDT = -(EUR_SPEED - QS); 
    
    EUR(:,:,1) = EUR(:,:,1) - Q1;
    EUR(:,:,2) = EUR(:,:,2) - Q2;
   
    
    EUR_Q_SPEED = sum(EUR.^2,3).^.5;
    EUR_Q_SPEED = imfilter(EUR_Q_SPEED,fspecial('disk',41),'replicate');
    NOR_EUR_S = bsxfun(@times,EUR,EUR_Q_SPEED.^-1);
    
    [gradS(:,:,2),gradS(:,:,1)] = gradient(EUR_Q_SPEEDT);
    STRAIN = sum(gradS.*NOR_EUR_S,3);
    %STRAIN = sum(gradS.^2,3).^.5;
    
    
    
    
    
    tmpS = ba_interp2(STRAIN,np1,np2);
    tmpV = ba_interp2(EUR_Q_SPEED,np1,np2);
    
    tmpZ = Z(:,:,e);
    tmpZ(idx) = tmpS;
    Z(:,:,e) = tmpZ;
      
    tmpZ1 = Z1(:,:,e);
    tmpZ1(idx) = tmpV;
    Z1(:,:,e) = tmpZ1;
    %{
    G1 = ba_interp2(gradS(:,:,1),np1,np2);
    G2 = ba_interp2(gradS(:,:,2),np1,np2);
    imshow(tmpZ1,[]);
    hold on
    for e = 1:100:numel(np1)
        plot([np1(e) np1(e)+1000*G2(e)],[np2(e) np2(e)+1000*G1(e)],'g')
        drawnow
    end
    %}
    %imshow(tmpZ,[]);drawnow
    e
end