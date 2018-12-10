function [] = stomata_cnnVersion(imageFile,convnet,freqK,R,oPath,rPath)
tic


if ~isempty(oPath)
    mkdir(oPath)
end
    
    
    
    

    [pth,nm,ext] = fileparts(imageFile);

    fft = extractSingleBugEye_v2(imageFile,freqK,R,'');
    
    
    
    tmp = permute(fft.f,[3 4 1 2]);
    tsz = size(tmp);
   % tmp = reshape(tmp,[tsz(1:2) 1 433*433]);
    tmp = reshape(tmp,[tsz(1:2) 1 473*473]);
   
    
    
    tmpA = fft.A;
    tmpA = diff(cat(4,tmpA,tmpA(:,:,:,end)),1,4);
    perE = abs(tmpA)/(2*pi);
    tmpA = .5*(((perE.^2)+1).^.5 + (((1-perE).^2)+1).^.5);
    
    
    tmpA = permute(tmpA,[3 4 1 2]);
    tsz2 = size(tmpA);
    %tmpA = reshape(tmpA,[tsz2(1:2) 1 433*433]);
    tmpA = reshape(tmpA,[tsz(1:2) 1 473*473]);
    
   
    
    DATA = cat(3,tmp,tmpA,zeros(size(tmp)));
    clear tmp tmp2
    [classLabel PROB] = classify(convnet,single(DATA));
    
   
    J = imread(imageFile);
    for r = 1:4
        J(1:20,:) = [];
        J = imrotate(J,90);
    end
    J = imresize(J,[473 473]);
    
    
    
    
    %classLabel = reshape(double(classLabel),[433 433]);
    classLabel = reshape(double(classLabel),[473 473]);
    
    
    mask = classLabel == 2;
    prob = reshape(double(PROB(:,2)),[473 473]);
    
    
    mask = prob >= .4;
    mask = imclose(mask,strel('disk',5,0));
    mask = bwareaopen(mask,50);
    mask = imfill(logical(mask),'holes');
    
    
    
    R = regionprops(mask==1,'Centroid','Area','Eccentricity','Perimeter','Orientation','MajorAxisLength');
    fidx = count([R.Area].^.5);
    fidx = count([R.MajorAxisLength]);
    fidx = find(fidx >= 2);
    
    % 
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
    
    
    
    
    image(cat(3,J,J,J))
    hold on
    axis off
    plot(cLOC(:,1),cLOC(:,2),'r*')
    plot(newLOC2(:,1),newLOC2(:,2),'g*')
    plot(bLOC(:,1),bLOC(:,2),'bo');
    
    %{
    out = flattenMaskOverlay(bindVec(J),logical(mask));
    imshow(out,[]);hold on
    for l = 1:numel(fidx)
       plot(R(fidx(l)).Centroid(1),R(fidx(l)).Centroid(2),'b*') 
    end
    %}
    
    %{
    
    prob1 = reshape(double(PROB(:,2)),[433 433]);
    prob1 = imfilter(prob1,fspecial('gaussian',[21 21],2),'replicate');
    prob1 = log(prob1);
    
    
    %prob1f = imfilter(prob1,fspecial('gaussian',[21 21],2),'replicate');
    %{
    %prob1f = -prob1f;
    for e = 1:100
        prob1f = prob1f + .3*del2(prob1f);
        %contour(prob1f,linspace(min(prob1f(:)),max(prob1f(:)),10));
        %drawnow
    end
    %prob1f = -prob1f;
     %}
    L = linspace(0,-10,100);
    %L = linspace(.02,1,100);
    cnt = 1;
    ECC = [];
    CEN = [];
    AREA = [];
    PER = [];
    for level = 1:numel(L)
        slice = prob1 > L(level);
        %imshow(slice);
        R = regionprops(slice,'Eccentricity','Centroid','Area','Perimeter');
        for p = 1:numel(R)
            CEN(cnt,:) = R(p).Centroid;
            ECC(cnt) = R(p).Eccentricity;
            PER(cnt) = R(p).Perimeter;
            AREA(cnt) = R(p).Area;
            cnt = cnt + 1;
        end
        %drawnow
    end
    AREA = (AREA/pi).^.5;
    rm = find(ECC > .7 | AREA < 4);
    %rm = [];
    ECC(rm) = [];
    AREA(rm) = [];
    CEN(rm,:) = [];
    PER(rm) = [];
    
    mag = 1;
    fidx = find(count2(ECC,mag).*count2(AREA,mag).*count2(PER,mag));
    cM = zeros(size(prob1));
    CEN = round(CEN);
    for c = 1:numel(fidx)
        cM(CEN(fidx(c),2),CEN(fidx(c),1)) = cM(CEN(fidx(c),2),CEN(fidx(c),1)) + 1;
    end
    
    cM = cM >= 1;
    %{
    cM = imfilter(cM,fspecial('gaussian',[21 21],5),'replicate');
    cM = (imdilate(cM,strel('disk',5,0)) == cM) & cM ~= 0;
    %}
    cM = imdilate(cM,strel('disk',3,0));
    %{
    dcM = imdilate(cM,strel('disk',5,0));
    
     
    
    
    prob1f = imfilter(prob1,fspecial('gaussian',[21 21],3),'replicate');
    for e = 1:100
        prob1f = prob1f + .5*del2(prob1f);
        contour(prob1f,linspace(min(prob1f(:)),max(prob1f(:)),10));
        drawnow
    end
    
    
    mask = imdilate(prob1f,strel('disk',5,0)) == prob1f;
    
    Iobrd = imdilate(prob1, strel('disk',7,0));
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    
    
   
    midx = find(mask);
    v = prob1(midx);
    rIDX = exp(v) < .1;
    mask(midx(rIDX)) = 0;
    %{
    
    %%%%%%%%%%%%%%%
    K = imfill(mask,'holes');
    K = imerode(K,strel('disk',3,0));
    K = logical(K);
    K = bwareaopen(K,50);
    K = myStomataPeaks(K);
    R = regionprops(K,'Centroid');
    Z = zeros(size(K));
    for i = 1:numel(R)
        Z(round(R(i).Centroid(2)),round(R(i).Centroid(1))) = 1;
    end
    Z = logical(imdilate(Z,strel('disk',7)));
    %%%%%%%%%%%%%%%
    %}
    %mask = mask .* (prob1f > .1);
    dmask = imdilate(mask,strel('disk',11,0));
    %mask = logical(classLabel==2);
    %mask = bwareaopen(mask,10);
   
    
    
    
    
    
    %}
    
    mask = cM;
    
    %}
   
    
    
    if ~isempty(oPath)
        fileList = {};
        fileList{end+1} = [oPath filesep nm '_labeledImage.jpg'];
        saveas(gca,fileList{end});

        fileList{end+1} = [oPath filesep nm '_UNlabeledImage.jpg'];
        imwrite(J,fileList{end});


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