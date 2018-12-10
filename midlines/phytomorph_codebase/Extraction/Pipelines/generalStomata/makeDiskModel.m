function [stomataProb] = makeDiskModel(toD,numPC)
    %% new SAMP disk
    [ssU,ssE,rrU,rrE] = dc(toD,numPC);



    %% T BUILD
    close all
    testC = [];
    dE = [];
    errD = [];
    for tr = 1:size(toD,3)
        [dE(:,tr),testC(:,tr),errD(:,tr)] = applyTC(toD(:,:,tr),ssU,ssE,rrU,rrE,[]);
    end

    clear modelStruct
    for d = 1:size(testC,1)
        [modelStruct(d).F,modelStruct(d).X] = ksdensity(testC(d,:));
        modelStruct(d).F = modelStruct(d).F / sum(modelStruct(d).F);
        figure;
        plot(modelStruct(d).X,modelStruct(d).F)
    end


    errD = reshape(errD,size(toD));
    [ssUe,ssEe,rrUe,rrEe] = dc(errD,[1 1]);
    for tr = 1:size(errD,3)
        [dEe(:,tr),testCe(:,tr),errDe(:,tr)] = applyTC(errD(:,:,tr),ssUe,ssEe,rrUe,rrEe,[]);
    end

    clear modelStructE
    for d = 1:size(testCe,1)
        [modelStructE(d).F,modelStructE(d).X] = ksdensity(testCe(d,:));
        modelStructE(d).F = modelStructE(d).F / sum(modelStructE(d).F);
        figure;
        plot(modelStructE(d).X,modelStructE(d).F)
    end

    U1{1} = ssU;
    E1{1} = ssE;
    U2{1} = rrU;
    E2{1} = rrE;
    U1{2} = ssUe;
    E1{2} = ssEe;
    U2{2} = rrUe;
    E2{2} = rrEe;
    MS{1} = modelStruct;
    MS{2} = modelStructE;
    stomataProb = @(X)recursiveTC(X,U1,E1,U2,E2,MS);
    %stomataProb = @(X)applyTC(X,U1,E1,U2,E2,modelStruct);
end

function [ssU,ssE,rrU,rrE] = dc(toD,numPC)
    close all
    zSAMP = [];
    %rSAMP = circshift(SAMP3,100,2);
    %pSAMP = cat(3,SAMP3,rSAMP);


    pSAMP = toD;
    sSZ = size(pSAMP);

    for e = 1:size(pSAMP,3)
        tmp = pSAMP(:,:,e);
        %tmp = zscore(tmp(:),1,1);
        tmp = reshape(tmp,sSZ(1:2));
        zSAMP(:,:,e) = tmp; 
    end

    % sim radial vectors
    nPC = 3;
    sSZ = size(zSAMP);
    zSAMP = reshape(zSAMP,[sSZ(1) prod(sSZ(2:3))]);
    [ssU,ssE,ssL] = PCA_FIT_FULL_Tws(zSAMP,numPC(1));
    zSAMP = reshape(zSAMP,sSZ);
    zSAMP = reshape(zSAMP,[sSZ(1) prod(sSZ(2:3))]);
    [ssC] = PCA_REPROJ_T(zSAMP,ssE,ssU);
    zSAMP = reshape(zSAMP,sSZ);
    ssC = reshape(ssC,[numPC(1) sSZ(2:3)]);


    plot(ssE);
    %waitforbuttonpress
    figure
    plot(mean(ssC,3)');
    %waitforbuttonpress

    % sim radial vectors
    nPCr = 3;
    zSAMP = permute(zSAMP,[2 1 3]);
    sSZ = size(zSAMP);
    zSAMP = reshape(zSAMP,[sSZ(1) prod(sSZ(2:3))]);
    [rrU,rrE,rrL] = PCA_FIT_FULL_Tws(zSAMP,numPC(2));
    zSAMP = reshape(zSAMP,sSZ);
    zSAMP = reshape(zSAMP,[sSZ(1) prod(sSZ(2:3))]);
    [rrC] = PCA_REPROJ_T(zSAMP,rrE,rrU);
    zSAMP = reshape(zSAMP,sSZ);
    rrC = reshape(rrC,[numPC(2) sSZ(2:3)]);
    zSAMP = ipermute(zSAMP,[2 1 3]);
    rrC = ipermute(rrC,[2 1 3]);
    figure;
    plot(rrE);
    %waitforbuttonpress
    figure
    plot(mean(rrC,3));
    %% make tensor
    sSZ = size(rrC);
    zTMP = reshape(rrC,[sSZ(1) prod(sSZ(2:3))]);
    [finalC] = PCA_REPROJ_T(zTMP,ssE,ssU);
    zTMP = reshape(zTMP,sSZ);
    finalC = reshape(finalC,[numPC(1) sSZ(2:3)]);
    szC = size(finalC);
    finalC = reshape(finalC,[prod(szC(1:2)) szC(3)]);
    [ccU,ccE,ccL] = PCA_FIT_FULL_Tws(finalC,size(finalC,1));
end