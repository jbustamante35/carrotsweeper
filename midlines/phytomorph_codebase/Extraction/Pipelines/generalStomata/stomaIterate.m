function [PARA] = stomaIterate(I,p,xD,trainedNetwork_Angle_local,Z0,trainedNetwork_uM,Z1,sz,funcGrade,disp)
    PARA = [0 1 1 0 0];

   
    
    
    p0 = p;
    [pCordinate(2),pCordinate(1)] = ind2sub(size(I),p);
    pCordinate = double(pCordinate);
    displacement = [0 0];
    rot = 0;
    loopMax = 10;
    w = @(l)(.5*exp((-1/(1.5*loopMax))*(l-1)));
    for loop = 1:20
        [T] = generateDomainTransformation(PARA);
        [xDs] = transformDomain(xD,T);
        [samp] = myInterp2Sampler(I,pCordinate',xDs,sz);
        sz = size(samp);
        samp = zscore(samp(:),1,1);
        samp = reshape(samp,sz);
        [paraDisplay,cP] = generatePARAfromNetworks(samp,trainedNetwork_Angle_local,Z0,trainedNetwork_uM,Z1);
        
        
        [paraUpdate,displacement,newPara] = suggestNewPara(paraDisplay,cP,sz,w(loop));
        
        if disp
            displayStoma(samp,paraDisplay,[31 31],cP','r');
            displayStoma([],newPara,[31 31],cP','b');
            displayStoma([],paraUpdate,[31 31],cP','g');
            drawnow
        end
        PARA(1) = PARA(1) + paraUpdate(1);
        waitforbuttonpress
        %lV(loop) = PARA(1);
        %rot = rot + paraUpdate(1);
        %title(num2str(loop));
        %waitforbuttonpress
        %drawnow
        %if loop == 1
            %waitforbuttonpress
        %end
    end


    


    if ~isempty(funcGrade)
        gPARA = double([PARA(1) paraUpdate(2:3) 0 0 ]);
        [g1,g2] = ndgrid(linspace(-1.5,1.5,50),linspace(-1.5,1.5,50));
        G = [g2(:),g1(:)];
        [T] = generateDomainTransformation(gPARA);
        [G] = transformDomain(G,T);
        [samp] = myInterp2Sampler(I,pCordinate',G,size(g1));
        sz = size(samp);
        samp = zscore(samp(:),1,1);
        samp = reshape(samp,sz);
        grade = funcGrade(samp);
        PARA = [PARA,grade(:)'];
    end

    


end