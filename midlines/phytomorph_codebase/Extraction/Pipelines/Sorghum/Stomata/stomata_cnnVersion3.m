function [] = stomata_cnnVersion2(imageFile,convnet,freqK,R,oPath,rPath,pSZ)
tic
HARDSIZE = 33;

if ~isempty(oPath)
    mkdir(oPath)
end
    
    
    
    diskSample = 1;

    if diskSample
        [R,T] = ndgrid(linspace(0,(pSZ-1)/2,(pSZ-1)/2),linspace(-pi,pi,round(2*pi*(pSZ-1)/2)));
        X1 = R.*cos(T) + (pSZ-1)/2;
        X2 = R.*sin(T) + (pSZ-1)/2;
        DISK = cat(3,X1,X2);
    end
    

    [pth,nm,ext] = fileparts(imageFile);
    I = imread(imageFile);
    SZ = size(I) - (HARDSIZE-1);
    sliding = false;
    if sliding
        B = colfilt(I,[HARDSIZE HARDSIZE],'sliding',@(block)lowMemCNN(block,convnet,[HARDSIZE HARDSIZE],2,DISK));
    else

        I = im2col(I,[HARDSIZE HARDSIZE],'sliding');
        I = reshape(I,[HARDSIZE HARDSIZE 1 size(I,2)]);
        
        diskSample = 1;
        if diskSample
            nI = zeros([size(DISK,1) size(DISK,2) size(I,3) size(I,4)]);
            parfor s = 1:size(nI,4)
                nI(:,:,:,s) = ba_interp2(I(:,:,:,s),X1,X2);
                %nI(:,:,:,s) = fft2(nI(:,:,:,s));
            end
            I = nI;
        end
    
        
        
        
        [classLabel PROB] = classify(convnet,I,'MiniBatchSize',2048*4);
        rPROB = reshape(PROB,[SZ(1) SZ(2) size(PROB,2)]);

        classLabel = reshape(double(classLabel),[SZ(1) SZ(2)]);


        J = imread(imageFile);
        for r = 1:4
            J(1:10,:) = [];
            J = imrotate(J,90);
        end
        J = imresize(J,[473 473]);
    end
    %B = bindVec(B);
    fB = imfilter(B,fspecial('disk',5));
    %fB = B;
    %marker = imerode(fB,strel('disk',3,0)) .* B >. 5;
    marker = imdilate(fB,strel('disk',7,0)) == fB;
    recon = imreconstruct(marker.*fB-.1,fB);
    mask = (recon == (fB-.1)) & (fB > .50);
    %recon = bindVec(recon);
    %mask = (imdilate(recon,strel('disk',7,0)) == recon) & (B > .6) & (recon > .48);
    
    %mask = bwareaopen(mask,10);
    
    
    R = regionprops(mask==1,'Centroid','Area','Eccentricity','Perimeter','Orientation','MajorAxisLength');
    fidx = count([R.Area].^.5);
    fidx = count([R.MajorAxisLength]);
    R(fidx == 0) = [];
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
    
    
    
 
    figure
    image(cat(3,I,I,I))
    hold on
    axis off
    plot(cLOC(:,1),cLOC(:,2),'r*')
    if ~isempty(newLOC2)
        plot(newLOC2(:,1),newLOC2(:,2),'g*')
        plot(bLOC(:,1),bLOC(:,2),'bo');
    end
    
   
    
    
    if ~isempty(oPath)
        fileList = {};
        fileList{end+1} = [oPath filesep nm '_labeledImage.jpg'];
        saveas(gca,fileList{end});

        fileList{end+1} = [oPath filesep nm '_UNlabeledImage.jpg'];
        imwrite(I,fileList{end});


        fileList{end+1} = [oPath filesep nm '_locations.csv'];
        csvwrite(fileList{end},cLOC);
        
        fileList{end+1} = [oPath filesep nm '_WITH2locations.csv'];
        csvwrite(fileList{end},[sLOC;newLOC2]);


        pushToiRods(rPath,fileList)
        
        close all
    end
    close all
    
    toc
end