function [MAJOR MINOR rowNumber KERD F ROT] = processCrossSection(I,disp)
    DS = 50;
    samplePER = 1.5;
    PER_SHRINK = .95;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make gray, mask and simple process mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = rgb2gray(I);
    M = G > graythresh(G);
    M = bwareaopen(M,10000);
    M = imclose(M,strel('disk',110));
    R = regionprops(M,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
    MAJOR = .5*R(1).MajorAxisLength*PER_SHRINK;
    MINOR = .5*R(1).MinorAxisLength*PER_SHRINK;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % crop the objects at PER greater than diameter 
    %OUTER = samplePER*max(sum(M,2))/2;
    OUTER = samplePER*R(1).MajorAxisLength;
    % make circular sample disk
    %[n1 n2] = ndgrid(linspace(-pi,pi,2000),linspace(0,OUTER,OUTER));
    ROT = [[cos(R(1).Orientation*pi/180) sin(R(1).Orientation*pi/180)];[-sin(R(1).Orientation*pi/180) cos(R(1).Orientation*pi/180)]];
    %{
    [n1 n2] = ndgrid(linspace(-pi,pi,2000),linspace(0,1,round(OUTER)));
    X = R(1).MajorAxisLength*n2.*cos(n1);
    Y = R(1).MinorAxisLength*n2.*sin(n1);
    GR = [ROT*[X(:) Y(:)]']';
    X = GR(:,1) + R(1).Centroid(1);
    Y = GR(:,2) + R(1).Centroid(2);
    %}
    
    NS = 3*round(.5*samplePER*R(1).MinorAxisLength);
    [X,Y,RAD,Xn,Yn,RADn] = newDrawE(.5*samplePER*R(1).MajorAxisLength,.5*samplePER*R(1).MinorAxisLength,2000,ROT,R(1).Centroid,NS);
    
    
    
    % if display - show the sample rings downsampled by DS
    if disp
        h1 = figure;
        imshow(I,[]);
        hold on
        for e = 1:DS:size(X,2)
            plot(X(:,e),Y(:,e),'b');
        end
        %{
        for e = 1:DS:size(Xn,2)
            plot(Xn(:,e),Yn(:,e),'r');
        end
        %}
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the rings which capture the kernels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    % sample gray in disk hood
    S = ba_interp2(G,Xn,Yn);
    % smooth and gradient
    smoothedS = imfilter(S,fspecial('disk',5),'replicate');
    [d1 d2] = gradient(smoothedS);
    sig = bindVec(sum(abs(d2)));
    thresh = graythresh(sig);
    sig = sig > thresh;
    fidx = find(sig);
    MAX_RED = fidx(end);
    MIN_RED = fidx(1);
    %}
    [~,MAX_RED] = min(abs(.5*R(1).MajorAxisLength*PER_SHRINK - RADn(:,1)));
    if disp
        imshow(I,[]);
        hold on
        plot(Xn(:,MAX_RED),Yn(:,MAX_RED),'r');
    end
    
    [Xk Yk RADk] = newDrawE(RADn(MAX_RED,1),RADn(MAX_RED,2),500,ROT,R(1).Centroid,NS);
    smoothedG = imfilter(G,fspecial('disk',5),'replicate');
    S = ba_interp2(smoothedG,Xk,Yk);
    %smoothedS = imfilter(S,fspecial('disk',21),'replicate');
    [d1 d2] = gradient(S);
    
    ratio = -log(abs(d1)).*log(abs(d2)).^-1;
    %ratio = log(abs(d2)).*log(abs(d1)).^-1;
    
    
    
    
    testU = mean(abs(ratio),1);
    testS = std(abs(ratio),1,1);
    testA = testU.*testS;
    testA = imfilter(testA,ones(1,50),'replicate');
    dA = gradient(testA);
    sdA = imfilter(dA,ones(1,100),'replicate');
    %sdA = dA;
    idx = find(sdA == imdilate(sdA,strel('disk',40)));
    rmidx = idx/size(testA,2) > .8;
    idx(rmidx) = [];
    %MIN = idx(end);
    [~,sidx] = max(sdA(idx));
    MIN = idx(sidx);
    %{
    rmidx = idx/size(testA,2) > .9 | idx/size(testA,2) < .5;
    idx(rmidx) = [];
    %}
    %MIN = idx(1);
    %{
    rmidx = idx/size(testA,2) > .9 | idx/size(testA,2) < .4;
    idx(rmidx) = [];
    MIN = idx(1);
    %}
    
    %{
    idx = find(testA == imdilate(testA,strel('disk',100)));
    rmidx = idx/size(testA,2) > .9 | idx/size(testA,2) < .4;
    idx(rmidx) = [];
    MIN = idx(1);
    %}
    
    
    sig = bindVec(sum(abs(d2)));
    thresh = graythresh(sig);
    sig = sig > thresh;
    fidx = find(sig);
    %MAX = fidx(end);
    %MIN = fidx(1);
    
    
    
    MAX = size(Xk,2);
    if disp
        imshow(I,[]);
        hold on
        plot(Xk(:,MIN),Yk(:,MIN),'r');
        plot(Xk(:,MAX),Yk(:,MAX),'y');
    end
    
    KERD = RADk(MAX,1) - RADk(MIN,1);
    %{
    if disp
        imshow(I,[]);
        hold on
        plot(X(:,MIN),Y(:,MIN),'r');
        plot(X(:,MAX),Y(:,MAX),'y');
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
figure;
imshow(S,[]);
figure;
imshow(abs(d1),[])
figure;
imshow(abs(d2),[])
figure;
plot(sum(S,1))
figure;
plot(sum(abs(d2)));
%}
%{
figure;
plot(sum(abs(d1)),'r');
hold on
plot(sum(abs(d2)),'g');
    
    
    [g1 g2] = gradient(F);
%}

    %{
    MAJOR_OUT = R(1).MajorAxisLength*mean(n2(:,MAX));
    MINOR_OUT = R(1).MinorAxisLength*mean(n2(:,MAX));
    MAJOR_IN = R(1).MajorAxisLength*mean(n2(:,MIN));
    MINOR_IN = R(1).MinorAxisLength*mean(n2(:,MIN));
    
    KERNEL_DEPTH_MAJOR = MAJOR_OUT - MAJOR_IN;
    KERNEL_DEPTH_MINOR = MINOR_OUT - MINOR_IN;
    
    KERNEL_DEPTH_AVERAGE = .5*(KERNEL_DEPTH_MAJOR+KERNEL_DEPTH_MINOR);
    
    MIN = MAX - round(KERNEL_DEPTH_AVERAGE);
    %}
    
    
    
    
    PER = .5;
    WID = .5*(MAX-MIN);
    CEN = round(.5*(MIN+MAX));
    WID = round(PER*WID);
    VEC = CEN + (-WID:WID);
    
    
    %F = S(:,VEC);
    
    
    
    
    %[g1 g2] = gradient(F);
    sG = imfilter(G,fspecial('disk',11),'replicate');
    S = ba_interp2(sG,Xk,Yk);
    %smoothedS = imfilter(S,fspecial('disk',11),'replicate');
    [~,g2] = gradient(S(:,MIN:MAX));
    gZ = g2 > 0;
    sig = (gZ.*g2);
    sig = bsxfun(@minus,sig,mean(sig,1));
    sig1 = mean(abs(fft(sig)),2);
    %sig1 = imfilter(sig1,ones(3,1),'replicate');
    gZ = g2 < 0;
    sig = (gZ.*g2);
    sig = bsxfun(@minus,sig,mean(sig,1));
    sig2 = mean(abs(fft(sig)),2);
    %sig2 = imfilter(sig2,ones(3,1),'replicate');
    F = mean([sig1,sig2],2);
    %sig(1:6) = 0;
    %{
    sig = mean(F,2);
    sig = sig - mean(sig);
    F = mean(abs(fft(sig,[],1)),2);
    %}
    F(1:6) = 0;
    %figure;
    %plot(F);
    [J,idx] = max(F(1:end/2));
    T = idx-1;
    rowNumber = T;
    if disp
        figure(h1);
        rowCount = T;
        title(num2str(rowCount));
        TH = linspace(-pi,pi,2000);
        SIG = 50*cos(TH*T);
        Ro = .5*(MIN+MAX);
        xl = (Ro+SIG).*cos(TH) + R(1).Centroid(1);
        yl = (Ro+SIG).*sin(TH) + R(1).Centroid(2);
        plot(xl,yl,'r','LineWidth',3);
    end
end