function [grade,particleList] = stomataMinObjectiveFunc(I,vD,vSZ,hD,hsz,xD,sz,graderFunc,P,oP,gradeVec,fixerPara,regressNetworkFunc,gradeRatio,disp,startPara,paraVec,imageGradFunc,statEnergyMix)
    

    % handle offset for simplex search
    paraVec = paraVec - fixerPara;
    currentPara = paraVec;
    % render the point
    pToSample = paraVec(4:5);
    paraVec(4:5) = 0;
    %paraVec(1) = -paraVec(1);


    viewPara = paraVec;
    vpatch = samplePatch(I,vD,pToSample+P,viewPara,vSZ);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % allow network to suggest para at the "new" location and the "new" para
    % non-rot version of torsion from geo model
    [networkParaList,networkPointList,dX,networkOutput,networkPatch] = regressNetworkFunc(pToSample+P);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    [dpatch0,rpatch0,IdE(1,1),IdE(1,2)] = sampleTransformedRotatedPatch(I,xD,hD,pToSample+P,paraVec,sz,hsz,graderFunc,gradeVec,imageGradFunc,gradeVec);
    %[dpatch1,rpatch1,IdE(2,1),IdE(2,2)] = sampleTransformedRotatedPatch(I,xD,hD,pToSample+P,networkParaList(1,:),sz,hsz,graderFunc,gradeVec,imageGradFunc,gradeVec);
    %[dpatch2,rpatch2,IdE(3,1),IdE(3,2)] = sampleTransformedRotatedPatch(I,xD,hD,pToSample+P,networkParaList(2,:),sz,hsz,graderFunc,gradeVec,imageGradFunc,gradeVec);
    dpatch1 = dpatch0;
    dpatch2 = dpatch0;
    rpatch1 = rpatch0;
    rpatch2 = rpatch0;
    IdE(2:3,1:2) = 0;


    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the patch for the particle-0
    % obtain the sample at the proposed location
    [Tx0] = generateDomainTransformation(paraVec);
    [xDs0] = transformDomain(xD,Tx0);
    [patch0] = myInterp2Sampler(I,(pToSample+P)',xDs0,sz);
    sz = size(patch0);
    patch0 = zscore(patch0(:),1,1);
    patch0 = reshape(patch0,sz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the patch for the particle-2
    % obtain the sample at the proposed location
    [Tx1] = generateDomainTransformation(networkParaList(1,:));
    [xDs1] = transformDomain(xD,Tx1);
    [patch1] = myInterp2Sampler(I,(pToSample+P)',xDs1,sz);
    sz = size(patch1);
    patch1 = zscore(patch1(:),1,1);
    patch1 = reshape(patch1,sz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the patch for the particle-2
    % obtain the sample at the proposed location
    [Tx2] = generateDomainTransformation(paraVec);
    [xDs2] = transformDomain(xD,Tx2);
    [patch2] = myInterp2Sampler(I,(pToSample+P)',xDs2,sz);
    sz = size(patch2);
    patch2 = zscore(patch2(:),1,1);
    patch2 = reshape(patch2,sz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obtain the sample at the proposed location  - contraint image
    hPara = paraVec;
    hPara(2:3) = 1;
    [Hx] = generateDomainTransformation(hPara);
    [hD] = transformDomain(hD,Hx);
    [Hpatch] = myInterp2Sampler(I,(pToSample+P)',hD,hsz);
    sz = size(Hpatch);
    Hpatch = zscore(Hpatch(:),1,1);
    Hpatch = reshape(Hpatch,sz);
    HmodelStatisticalEnergy = imageGradFunc(Hpatch);
    HmodelStatisticalEnergy(isinf(-log(HmodelStatisticalEnergy))) = eps;
    HmodelStatisticalEnergy = sum(-log(HmodelStatisticalEnergy([1:2])));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grade the image on the model - likihood
    modelStatisticalEnergy = graderFunc(patch0(:)');
    modelStatisticalEnergy(isinf(-log(modelStatisticalEnergy))) = eps;
    modelStatisticalEnergy = sum(-log(modelStatisticalEnergy(gradeVec)));
    if modelStatisticalEnergy < eps
        fprintf('Warning: Energy @ zero K.\n');
        modelStatisticalEnergy = eps;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    StatisticalEnergyVector = IdE;
    integrationVec = size(StatisticalEnergyVector,1)^-1*ones(1,size(StatisticalEnergyVector,1));
    integrationVec = [1 0 0];
    %StatisticalEnergyVector = [modelStatisticalEnergy HmodelStatisticalEnergy];
    TotalStatisticalEnergy = integrationVec*StatisticalEnergyVector*statEnergyMix;
    %[TotalStatisticalEnergy,midx] = min(StatisticalEnergyVector*statEnergyMix);
    midx = 1;
    %{
    % agree here or somewher else nearby
    % travel cost too
    [dE bD] = bendEnergy(currentPara,startPara,pToSample+P-oP);
    [dE bD] = bendEnergy(currentPara,paraG,pToSample+P-oP);
    %}
    % energy between network paras
    [dE(1),bD(1,:)] = bendEnergy(networkParaList(1,:),networkParaList(2,:),currentPara(4:5)-startPara(4:5));
    [dE(2),bD(2,:)] = bendEnergy(networkParaList(1,:),currentPara,currentPara(4:5)-startPara(4:5));
    [dE(3),bD(3,:)] = bendEnergy(networkParaList(2,:),currentPara,currentPara(4:5)-startPara(4:5));
    totalEnergy(1,:) = [TotalStatisticalEnergy dE(1) TotalStatisticalEnergy+dE(1)];
    totalEnergy(2,:) = [TotalStatisticalEnergy dE(2) TotalStatisticalEnergy+dE(2)];
    totalEnergy(3,:) = [TotalStatisticalEnergy dE(3) TotalStatisticalEnergy+dE(3)];
    integrationVec = size(totalEnergy,1)^-1*ones(1,size(totalEnergy,1));
    deltaVarR = [1 0];
    grade = deltaVarR(1)*integrationVec*totalEnergy*gradeRatio + deltaVarR(2)*std(dE)^2;

    % return the particles that are making the choice
    networkParaList(1,4:5) = currentPara(4:5);
    networkParaList(2,4:5) = currentPara(4:5);
    particleList = [currentPara;networkParaList];
    particleList = [particleList;particleList(midx,:)];
    fprintf(['*********************************************\n'])
    fprintf(['*********************************************\n'])
    fprintf(['Total Energy Matrix:\n']);
    for e = 1:size(totalEnergy,1)
        fprintf([ num2str(totalEnergy(e,:),4) '\n'])
    end
    fprintf(['---------------------------------------------\n'])
    fprintf(['Energy Break Down:\n']);
    for e = 1:size(bD,1)
        fprintf([ num2str(bD(e,:),4) '\n'])
    end
    fprintf(['---------------------------------------------\n'])
    fprintf(['Statistical Energy Break Down:\n']);
    for e = 1:size(StatisticalEnergyVector,1)
        fprintf([ num2str(StatisticalEnergyVector(e,:),4) '\n'])
    end
    fprintf(['---------------------------------------------\n'])
    fprintf(['Grand Energy is:\n']);
    fprintf([num2str(grade) '\n']);
    fprintf(['---------------------------------------------\n'])
    fprintf(['Location:\n']);
    fprintf([num2str(paraVec(1:3)) '\n']);
    fprintf(['*********************************************\n'])
    fprintf(['*********************************************\n'])
    
    if disp
        CL = {'r' 'b' 'c'};
        boxP1 = [1 1 61 61];
        boxP2 = [61*3*1 1 61 61];
        boxP1(1:2) = boxP1(1:2) + [(midx-1)*61 0];
        boxP2(1:2) = boxP2(1:2) + [(midx-1)*61 0];
        
        oPara = startPara;
        oPara(2:3) = 1;
        oPoint = oPara(4:5);
        oPara(4:5) = 0;
        [O1,O2] = ndgrid(linspace(-30,30,61),linspace(-30,30,61));
        oD = [O2(:),O1(:)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obtain the sample at the proposed location
        [oT] = generateDomainTransformation(oPara);
        [oD] = transformDomain(oD,oT);
        [Opatch] = myInterp2Sampler(I,oPoint',oD,size(O1));
        sz = size(Opatch);
        Opatch = zscore(Opatch(:),1,1);
        Opatch = reshape(Opatch,sz);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%
        deformedParticlePatch = cat(1,dpatch0',dpatch1',dpatch2')';
        rotatedParticlePatch = cat(1,rpatch0',rpatch1',rpatch2')';
        deformedParticlePatch = imresize(deformedParticlePatch,size(rotatedParticlePatch));

        [d1,d2] = ndgrid(linspace(-2,2,size(rpatch0,1)),linspace(-2,2,size(rpatch0,2)));
        dsk = (d1.^2 + d2.^2).^.5;
        dsk = dsk < 1;
        dsk = bwboundaries(dsk);

        sz1 = size(deformedParticlePatch);
        sz2 = size(rotatedParticlePatch);
        sz3 = size(Opatch);
        sz4 = size(vpatch);



        I = cat(1,zeros(size(deformedParticlePatch,1),size(I,2)),I);
        %sz4 = size(Opatch);
        I(1:sz1(1),1:sz1(2)) = bindVec(deformedParticlePatch);
        I(1:sz2(1),(sz1(2)+1):(sz1(2)+sz2(2))) = bindVec(rotatedParticlePatch);
        I(1:sz3(1),(sz1(2)+sz2(2)+1):(sz1(2)+sz2(2)+sz3(2))) = bindVec(Opatch);
        I(1:sz4(1),(sz1(2)+sz2(2)+sz3(2)+1):(sz1(2)+sz2(2)+sz3(2)+sz4(2))) = bindVec(vpatch);

        P(2) = P(2) + size(deformedParticlePatch,1);
        startPara(end) = startPara(end) + size(deformedParticlePatch,1);
        paraDisplay = paraVec';
        displayStoma(I,paraDisplay,pToSample+P,'','r','.',[1 1]);


        toViewPoint = startPara(4:5);
        startPara(4:5) = 0;
        displayStoma('',startPara,toViewPoint,'','g','.',[1 1]);


        paraDisplay = networkParaList(1,:);
        displayStoma('',paraDisplay,pToSample+P,'','b','.',[1 1]);

        paraDisplay = networkParaList(2,:);
        displayStoma('',paraDisplay,pToSample+P,'','c','.',[1 1]);

        %plot(dsk{1}(:,1),dsk{1}(:,2),'r')

        rectangle('Position',boxP1,'EdgeColor',CL{midx});
        rectangle('Position',boxP2,'EdgeColor',CL{midx});
        plot([1 61*6],[31 31],'m')
        drawnow
        hold off
    end

end

function [dE dEBreakDown] = bendEnergy(para1,para2,displacement)
    twist =3* norm(para1(1) - para2(1));
    stretchMajor = norm(para1(2) - para2(2));
    stretchMinor = norm(para1(3) - para2(3));
    twist = 0;
    stretchMajor = 0;
    stretchMinor = 0;
    %displacement = norm(para1(4:5) - para2(4:5));
    displacement = norm(displacement);
    %displacement = 0;
    eVec = [twist stretchMajor stretchMinor displacement];
    
    dE = sum(eVec);
    dEBreakDown = eVec;
end


function [patch] = samplePatch(I,domain,point,tranformPara,sz)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the patch for the particle-2
    % obtain the sample at the proposed location
    [T] = generateDomainTransformation(tranformPara);
    [Ddomain] = transformDomain(domain,T);
    [patch] = myInterp2Sampler(I,point',Ddomain,sz);
    sz = size(patch);
    patch = zscore(patch(:),1,1);
    patch = reshape(patch,sz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [dpatch,rpatch,DdE,RdE] = sampleTransformedRotatedPatch(I,Ddomain,Rdomain,point,tranformPara,szD,szR,gradeFuncD,Didx,gradeFuncR,Ridx)
    dpatch = samplePatch(I,Ddomain,point,tranformPara,szD);
    tranformPara(2:3) = 1;
    rpatch = samplePatch(I,Rdomain,point,tranformPara,szR);
    DdE = gradePatch(dpatch,gradeFuncD,Didx);
    RdE = gradePatch(rpatch,gradeFuncR,Ridx);
end


function [dE] = gradePatch(patch,gradeFunc,gradeIDX)
    dE = gradeFunc(patch);
    dE(isinf(-log(dE))) = eps;
    dE = -log(dE(gradeIDX))
    dE = sum(dE);
    %dE = sum(-log(dE));
    if dE < eps
        fprintf('Warning: Energy @ zero K.\n');
        modelStatisticalEnergy = eps;
    end
end
