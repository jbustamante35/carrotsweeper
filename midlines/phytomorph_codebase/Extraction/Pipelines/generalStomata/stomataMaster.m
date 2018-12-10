function [] = stomataMaster(I,trainedNet,R,Domain,N,freqToKeep,newE,newU,newSZ,oPath,rPath)




    plFunc = @(X)generateImageDomain(X,R(2));
    % generate the sampler function
    samplerFunction = @(X,P)myInterp2Sampler(X,P,Domain,N);
    % generate the preprocess function
    preprocessFunction = @(X,P)fftPatch(X,P,samplerFunction,freqToKeep);
    % generate the network apply

    networkApply = @(X,P)applyNetworkToPoint(X,P,preprocessFunction,trainedNet,newE,newU);



    if ischar(I)
        [pth nm ext] = fileparts(I);
        I = double(imread(I))/255;
    end
    
 
    
    networkResults = applyFuncToLocation(I,networkApply,plFunc,0,1);
    networkResults = squeeze(networkResults);
    networkResults = reshape(networkResults,[newSZ size(networkResults,2)]);
    
    
    toRemove = (size(I) - newSZ)/2;
    J = I;
    for e = 1:4
        J(:,1:toRemove) = [];
        J = imrotate(J,90);
    end
    
    [value,mask] = max(networkResults,[],3);
    mask = mask == 3 & value > .75;
    
    %mask = networkResults(:,:,3) > .5;
    
    mask = bwareaopen(mask,60);
    %mask = imclose(mask,strel('disk',5,0));
    %mask = bwareaopen(mask,50);
    mask = imfill(logical(mask),'holes');
    
    
    
    R = regionprops(mask==1,'Centroid','Area','Eccentricity','Perimeter','Orientation','MajorAxisLength');
    fidx = count([R.Area].^.5);
    fidx = count([R.MajorAxisLength]);
    fidx = find(fidx >= 2);
    
    % 
    cLOC = [];
    for e = 1:numel(R)
        cLOC(e,:) = R(e).Centroid;
    end
    cLOC2 = cLOC(fidx,:);
    sLOC = cLOC;
    bLOC = cLOC(fidx,:);
    sLOC(fidx,:) = [];
    theta = [R(fidx).Orientation];
    L = [R(fidx).MajorAxisLength]/2;
    per = .75;
    dX = per*L.*cos(theta*pi/180);
    dY = -per*L.*sin(theta*pi/180);
    newLOC2 = [];
    for e = 1:size(cLOC2,1)
        newLOC2 = [newLOC2 ; [cLOC2(e,:) + [dX(e) dY(e)]] ; [cLOC2(e,:) - [dX(e) dY(e)]]];
    end
    
    
    close all
    out = flattenMaskOverlay(255*J,mask);
    image(out)
    hold on
    axis off
    plot(cLOC(:,1),cLOC(:,2),'r*')
    if ~isempty(newLOC2)
        plot(newLOC2(:,1),newLOC2(:,2),'g*')
        plot(bLOC(:,1),bLOC(:,2),'bo');
    end
    
    
     if ~isempty(oPath)
         
         
        mkdir(oPath);
        fileList = {};
        fileList{end+1} = [oPath filesep nm '_labeledImage.jpg'];
        saveas(gca,fileList{end});

        fileList{end+1} = [oPath filesep nm '_UNlabeledImage.jpg'];
        imwrite(J,fileList{end});


        fileList{end+1} = [oPath filesep nm '_locations.csv'];
        csvwrite(fileList{end},cLOC);
        
        fileList{end+1} = [oPath filesep nm '_WITH2locations.csv'];
        csvwrite(fileList{end},[sLOC;newLOC2]);

        fileList{end+1} = [oPath filesep nm '_probs.tif'];
        imwrite(networkResults,fileList{end});

        pushToiRods(rPath,fileList)
        
        close all
    end
    close all
    
    
end