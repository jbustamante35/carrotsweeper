function [circleMatrix MASK tMASK OFFSET] = getCircles(imageStack,rec,N,topN)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load, grayscale, and transform N images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CS = [];
    toColorStack = 1;
    disp = 0;
    parfor e = 2:N
        fprintf(['start: reading, transforming image:' num2str(e) '\n']); 
        tmpI = double(imread(imageStack{e}))/255;
        if toColorStack
            CSS(:,:,:,e-1) = tmpI;
        end
        %tmpI = rgb2gray(tmpI);
        tmpI = rgb2hsv_fast(tmpI,'','H');
        tmpI = tmpI(:,:,1);
        oSZ = size(tmpI);
        [tmpI dx dy] = imtransform(tmpI,rec,'FillValues',1);
        tSZ = size(tmpI);
        CS(:,:,e-1) = tmpI;
        fprintf(['end: reading, transforming image:' num2str(e) '\n']);
    end
    uCS = mean(CS,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load, grayscale, and transform N images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transform one image to get transformation values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpI = imread(imageStack{1});
    %tmpI = rgb2gray(tmpI);
    oSZ = size(tmpI);
    [tmpI dx dy] = imtransform(tmpI,rec,'FillValues',0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transform one image to get transformation values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obtain mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get average image
    %{
    M = double(uCS)/255 > .8;
    M = bwareaopen(M,50);
    M = bwmorph(M,'skel','inf');
    M = imdilate(M,strel('disk',4,0));
    OZ = round([340 385]/2);
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obtain mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hough transform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: hough transform \n']);
    OZ = round([165 190]/2);
    OZ = round([125 160]/2);
    M = uCS < .06;
    M = bwareaopen(M,200);
    copyVec = M(end,:);
    M(end,:) = 0;
    M = imclearborder(M);
    M(end,:) = copyVec;
    M = imfill(M,'holes');
    [centers, radii, metric] = imfindcircles(M,OZ,'ObjectPolarity','bright','Sensitivity',.99);
    mini = 4;
    CM = [];
    RM = [];
    for e = 1:size(centers,1)
        try
            CM(round(centers(e,2)/4),round(centers(e,1)/4)) = metric(e) + rand(1)*max(metric)*.0001;;
            RM(round(centers(e,2)/4),round(centers(e,1)/4)) = radii(e);
        catch
            
        end
    end
    fprintf(['ending: hough transform \n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hough transform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % non max suppress circles based on hough strength
    % select out the top 168 circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: non-max suppression \n']);
    miniCM = imdilate(CM,strel('disk',round(1.75*mean(radii)/mini),0));
    cp = (miniCM == CM).*(CM ~= 0);
    cp_index = find(cp);
    cp_strength = CM(cp_index);
    [J sidx] = sort(cp_strength,'descend');
    [cp1 cp2] = find(cp);
    %{
    if numel(cp1) < 100
        topN = topN / 2;
    end
    %}
    cp_index = cp_index(sidx(1:topN));
    cp1 = cp1(sidx(1:topN));
    cp2 = cp2(sidx(1:topN));
    fprintf(['ending: non-max suppression \n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % non max suppress circles based on hough strength
    % select out the top 185 circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if disp
        imshow(uCS,[]);
        viscircles(mini*[cp2 cp1],  1.2*RM(cp_index),'EdgeColor','r');
        drawnow
    end
    
    rad_vec = 1.2*RM(cp_index);
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the average gray scale value within each circle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: disk averaging \n']);
    rad_vec = 1.2*RM(cp_index);
    circleMean = [];
    parfor e = 1:numel(cp1)
        fprintf(['starting: disk average:' num2str(e) '\n']);
        ZM = zeros(size(M));
        ZM(mini*cp1(e),mini*cp2(e)) = 1;
        ZM = bwdist(ZM);
        ZM = ZM < rad_vec(e);
        circleMean(e) = mean(double(uCS(find(ZM))));
        fprintf(['ending: disk average:' num2str(e) '\n']);
    end
    fprintf(['starting: disk averaging \n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the average gray scale value within each circle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the brightest circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    topN = 168;
    [J sidx] = sort(circleMean,'descend');
    cp1 = cp1(sidx(1:topN));
    cp2 = cp2(sidx(1:topN));
    cp_index = cp_index(sidx(1:topN));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the brightest circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    if disp
        imshow(uCS,[]);
        viscircles(mini*[cp2 cp1],  1.2*RM(cp_index),'EdgeColor','r');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the master mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ZM = zeros(size(M));
    for e = 1:numel(cp1);
        ZM(mini*cp1(e),mini*cp2(e)) = 1;
    end
    ZM = bwdist(ZM);
    MASK = ZM < mean(rad_vec)/1.2;
    tMASK = MASK;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the master mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the master mask inv transform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [X Y] = find(MASK);
    X = X + dy(1);
    Y = Y + dx(1);
    OFFSET = [dx(1) dy(1)];
    [U,V] = tforminv(rec,Y,X);
    MASK = zeros(oSZ(1:2));
    for e = 1:numel(U)
        MASK(max(round(V(e)),1),max(round(U(e)),1)) = 1;
    end
    MASK = imclose(MASK,strel('disk',5,0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the master mask inv transform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    
    circleMatrix = [mini*[cp2 cp1] 1.2*RM(cp_index)];
    
    
end