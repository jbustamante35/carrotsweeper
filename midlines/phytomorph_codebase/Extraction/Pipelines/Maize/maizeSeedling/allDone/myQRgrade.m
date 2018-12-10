function [grade,parts,value,RP] = myQRgrade(deltaPT,offsetPT,I,baseBOX,pt,NETS,offset,dpara)
    try
        deltaPT;
        newOFF = pt - offsetPT;
        pt = bsxfun(@plus,deltaPT,newOFF);

        rm = (pt(:,1) < .5);
        pt(rm,1) = .5;
        
        %pt = pt + deltaPT - offsetPT;
        grade = 1;
        pt;

        MOD_baseBOX = (baseBOX'*pt(:,1)')';




        %[BOX] = point2Box(pt(:,2:3),repmat(MOD_baseBOX,[size(pt,1) 1]));
        for e = 1:size(MOD_baseBOX,1)
            [BOX(e,:)] = point2Box(pt(e,2:3),MOD_baseBOX(e,:));
        end
        BOX = round(BOX);


      
        
        
        
        EP = (BOX(:,1) + BOX(:,3)) ;
        rm = find(EP > size(I,2));
        BOX(rm,3) = BOX(rm,3) + size(I,2) - EP;
        
        EP = (BOX(:,2) + BOX(:,4));
        rm = find(EP > size(I,1));
        BOX(rm,4) = BOX(rm,4) + size(I,1) - EP;


        rm = find(BOX(:,1) <= 0);
        BOX(rm,1) = 1;
        rm = find(BOX(:,2) <= 0);
        BOX(rm,2) = 1;
        
        rm = find(BOX(:,3) <= 0);
        BOX(rm,3) = 100;
        rm = find(BOX(:,4) <= 0);
        BOX(rm,4) = 100;
        
        
        

        [FeatureSet,subI] = myCropper(I,BOX,.1,false,baseBOX);


        for r = 1:size(subI,4)
            subIN(:,:,:,r) = imresize(subI(:,:,:,r),[67 83]);

            subIN2(:,:,:,r) = imresize(imresize(subI(:,:,:,r),fliplr(baseBOX+1)),.08);
        end

        subI = subIN;
        %{
        subI = imcrop(I,BOX);

        %}

        boxCross = norm(BOX(3:4));
        RP = NETS{3}.predict(subIN2);
        RP = reshape(RP',[8 2 size(RP,1)]);
        cross1 = norm(RP(1,:) - RP(8,:));
        cross2 = norm(RP(3,:) - RP(6,:));
        cross = .5*(cross1 + cross2);
        boxCross = abs(cross - boxCross);
        boxCross = normpdf(boxCross,100,1000);
        boxCross0 = normpdf(100,100,1000);






        FeatureSet = cat(3,FeatureSet{1},FeatureSet{2});
        [type,prob] = NETS{1}.classify(FeatureSet);


        size(subI);


        scale = NETS{2}.predict(subI);

        scale = normpdf(scale,0,dpara(2)^2);
        scaleO = normpdf(dpara(1),0,dpara(2)^2);


        parts(:,1) = -scale.*scaleO.^-1;
        parts(:,2) = -prob*(1:3)';
        parts(:,3) = -boxCross.*boxCross0.^-1;
        value(:,1) = scale + offset;
        value(:,2:4) = prob;

        grade = sum(parts,2);
        %grade = parts(2);
        %grade = -(scale/scaleO + prob*(1:3)');
    catch ME
        ME
    end
    
    
end