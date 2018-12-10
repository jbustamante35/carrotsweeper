function [msg,qrCropBox] = getQRcode(I)
    %
    msg = '';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % threshold on Hue
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CL = reshape(I,[size(I,1)*size(I,2) 3]);
    hardMean = [0.8249    0.4369    0.3684];
    hardCov = [[0.0085    0.0094    0.0087];[0.0094    0.0109    0.0101];[0.0087    0.0101    0.0095]];
    P = mvnpdf(CL,hardMean,hardCov);
    P = reshape(P,[size(I,1) size(I,2)]);
    M = log(P) > -20;
    %{
    H = rgb2hsv_fast(I,'single','H');
    hsvI = rgb2hsv(I);
    H = hsvI(:,:,1);
    % also changed the threshold value - OLD .04
    val = .044;
    M = H < val | H > (1 - val);
    %}
    fM = imfill(M,'holes');
    fM = bwareaopen(fM,2000);
    fM = imdilate(fM,ones(50));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the largest object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = regionprops(fM,'Area','BoundingBox');
    [J,midx] = max([R.Area]);
    if isempty(R(midx).BoundingBox)
        return
    end
    dataCube = imcrop(I,R(midx).BoundingBox);
    qrCropBox = R(midx).BoundingBox;
    dataCubeMask = imcrop(M,R(midx).BoundingBox);
    
    % added to fill in the gaps from the bad threshold values
    % this was added at the phenome 2018 conference
    dataCubeMask = bwareaopen(dataCubeMask,200);
    % added to handle text in red box - after phenome 18
    dataCubeMask = imclose(dataCubeMask,strel('disk',31,0));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the third largest object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataCubeMaskf = imfill(dataCubeMask,'holes') - dataCubeMask;
    Rc = regionprops(logical(dataCubeMaskf),'Area','BoundingBox');
    [J,midxC] = sort([Rc.Area],'descend');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the qr data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qrCube = imcrop(dataCube,Rc(midxC(3)).BoundingBox);
    qrCube = rgb2gray(qrCube);
    %{
    for e = 1:4
        qrCube(1:10,:) = [];
        qrCube = imrotate(qrCube,90);
    end
    %}
    qrCube = imadjust(qrCube);
    qrCube = imsharpen(qrCube,'Amount',2);
    rotValue = linspace(0,360,360);
    
    
    for e = 1:numel(rotValue)
        qrCube_read = imrotate(qrCube,rotValue(e));
        try
            msg = decode_qr(qrCube_read);
        catch
            msg = [];
        end
        if ~isempty(msg)
            break
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the day cube image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dayCube = imcrop(dataCube,Rc(midxC(1)).BoundingBox);
    V = rgb2hsv_fast(dayCube,'single','V');
    H = rgb2hsv_fast(dayCube,'single','H');
    % change to .50 from .6 for carrot
    squareMask = H > .45 & V < .78;
    squareMask = bwareaopen(squareMask,100);
    squareMask = imclose(squareMask,strel('disk',5));
    squareMask = imfill(squareMask,'holes');
    BOXMask = edge(squareMask);
    squareMask = imclearborder(squareMask);
    blueMask = H > .53 & H < .91;
    blueMask = bwareaopen(blueMask,400);
    blueMask = imclose(blueMask,strel('disk',5));
    %{
    
    %blueMask = imdilate(blueMask,strel('disk',3));
    
    blueMaskS = bwmorph(blueMask,'skeleton',inf);
    blueMaskS = imdilate(blueMaskS,strel('disk',2));
    % find vertical linesV
    [H, theta, rho] = hough(blueMaskS,'Theta',linspace(-10,10,20));
    P  = houghpeaks(H,12,'Threshold',0);
    linesV = houghlines(blueMaskS,theta,rho,P,'FillGap',200,'MinLength',250);
    % find horizontal linesH
    [H, theta, rho] = hough(blueMaskS','Theta',linspace(-10,10,20));
    P  = houghpeaks(H,8,'Threshold',0);
    linesH = houghlines(blueMaskS',theta,rho,P,'FillGap',200,'MinLength',600);
    
    %{
    imshow(dayCube,[]);
    hold on
    for k = 1:length(linesV)
       xy = [linesV(k).point1; linesV(k).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
       % Plot beginnings and ends of linesV
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    end
    for k = 1:length(linesH)
       xy = [fliplr(linesH(k).point1); fliplr(linesH(k).point2)];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
       % Plot beginnings and ends of linesV
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    end
    %}
    
    BOXMask = zeros(size(blueMask));
    for k = 1:length(linesV)
       xy = [linesV(k).point1; linesV(k).point2];
       X = round(linspace(xy(1,1),xy(2,1),3000));
       Y = round(linspace(xy(1,2),xy(2,2),3000));
       for l = 1:numel(X)
           BOXMask(Y(l),X(l)) = 1;
       end
    end
    for k = 1:length(linesH)
       xy = [fliplr(linesH(k).point1); fliplr(linesH(k).point2)];
       X = round(linspace(xy(1,1),xy(2,1),3000));
       Y = round(linspace(xy(1,2),xy(2,2),3000));
       for l = 1:numel(X)
           BOXMask(Y(l),X(l)) = 1;
       end
    end
    
    %}
    BOXMask = imdilate(BOXMask,strel('disk',7));
    blueMask = blueMask.*BOXMask;
    
    blueMask = imdilate(blueMask,strel('disk',3));
    markerStrip = V < .6 - blueMask;
    markerStrip = bwareaopen(markerStrip,100);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the centers of each day cube
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ERO = 50;
    S1 = std(squareMask,1,1);
    bk1 = imerode(S1,ones(1,ERO));
    bk1 = imfilter(bk1,fspecial('disk',51));
    S1 = imfilter(S1,fspecial('disk',31));
    [lm1] = nonmaxsuppts(S1,81);
    R1 = regionprops(lm1,'Centroid');
    lm1 = zeros(size(lm1));
    for e = 1:numel(R1)
        lm1(round(R1(e).Centroid(1))) = 1;
    end
    S1 = bindVec(S1);
    level = graythresh(S1);
    lm1 = lm1 & S1 > level;
    
    S2 = std(squareMask,1,2);
    bk2 = imerode(S2,ones(ERO,1));
    bk2 = imfilter(bk2,fspecial('disk',51));
    S2 = imfilter(S2,fspecial('disk',31));
    [lm2] = nonmaxsuppts(S2,81);
    R2 = regionprops(lm2,'Centroid');
    lm2 = zeros(size(lm2));
    for e = 1:numel(R2)
        lm2(round(R2(e).Centroid(2))) = 1;
    end
    S2 = bindVec(S2);
    level = graythresh(S2);
    lm2 = lm2 & S2 > level;
   
    M = double(lm2)*double(lm1);
    [c1 c2] = find(M);
    centers = [c1 c2];
    %{
    imshow(dayCube,[]);
    hold on
    plot(c2,c1,'*');
    for e = 1:numel(c1)
        plot(c2(e),c1(e),'*');
        text(c2(e),c1(e),num2str(e));
    end
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % snap each center to nearest centroid 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = regionprops(squareMask,'Centroid','BoundingBox');
    CEN = [];
    for e = 1:numel(R)
        CEN = [CEN;fliplr(R(e).Centroid)];
    end
    for e = 1:size(centers,1)
        idx = snapTo(CEN,centers(e,:));
        boxTmp = imcrop(squareMask,R(idx).BoundingBox);
        markerTmp = imcrop(markerStrip,R(idx).BoundingBox);
        checked(e) = sum(boxTmp(:)-markerTmp(:))/sum(boxTmp(:)) < .98;
    end
    dayMatrix = [4:27];
    didx = find(checked);
    pictureDay = dayMatrix(didx(end));
    msg = [';' msg ';PictureDay:' num2str(pictureDay) ';'];
    
    msg = strrep(msg,';',';;');
    fidx = strfind(msg,';');
    msg(fidx(1:2:end)) = '}';
    msg(fidx(2:2:end)) = '{';
    msg(1) = [];
    msg(end) = [];
    msg = strrep(msg,':','_');
    msg = strrep(msg,'/','-');
    
end