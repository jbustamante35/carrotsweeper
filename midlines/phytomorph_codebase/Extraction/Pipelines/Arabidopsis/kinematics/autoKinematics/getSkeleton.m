function [skeleton] = getSkeleton(I,disp)
    SCALE = .25;
    I = imfilter(double(I),fspecial('disk',21),'replicate');
    tmpI = double(I)/255;
    [p1 p2] = gradient(tmpI);
    tmpI = (p1.^2 + p2.^2).^.5;
    % 500 200
    
    R = round([0 200]*SCALE);
    %R = round([0 60]*SCALE);
    N = round([100 R(2)*2*pi]*SCALE);
    tmpO = tmpI;
    tmpI = imresize(tmpI,SCALE);
    tmpR = tmpI;
    for r = 1:4
        tmpR(:,1:R(2)) = [];
        tmpR = imrotate(tmpR,90);
    end

    [n1 n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
    [d1 d2] = ndgrid(R(2):(size(tmpI,1)-R(2)),R(2):(size(tmpI,2)-R(2)));
    X = n1.*cos(n2);
    Y = n1.*sin(n2);
    
    NF = 20;
    g1 = {};
    g2 = {};
    g3 = {};
    fprintf(['**********************************\n']);
    fprintf(['Starting FFT:']);
    parfor p = 1:numel(d1)
        LOC = [d1(p) d2(p)];
        Xp = X + LOC(2);
        Yp = Y + LOC(1);
        if disp
            imshow(tmpI);
            hold on
            for s = 1:1:size(Xp,2)
                plot(Xp(:,s),Yp(:,s),'r');
            end
            for s = 1:1:size(Xp,1)
                plot(Xp(s,:),Yp(s,:),'r');
            end
            hold off
            drawnow
        end
        F = ba_interp2(tmpI,Xp,Yp);
        F = bsxfun(@minus,F,mean(F,2));
        fF = fft(F,[],2);
        ufF = mean(abs(fF),1);
        [q1 q2] = sort(ufF(1:round(end/2)),'descend');
        g1{p} = q1(1:NF);
        g2{p} = q2(1:NF);
        g3{p} = ufF(1:NF);
        fprintf(['.']);
        if mod(p,300) == 0
            fprintf('\n');
        end
    end
    fprintf('\n');
    fprintf(['Ending FFT:\n']);
    fprintf(['**********************************']);
    g1 = cell2mat(g1');
    g2 = cell2mat(g2');
    g3 = cell2mat(g3');
    T = reshape(g2(:,1:NF),[size(d1) NF]);
    AMP = reshape(g3(:,1:NF),[size(d1) NF]);
    sAMP = reshape(g1(:,1:NF),[size(d1) NF]);
    
    
    ampMASK = sAMP(:,:,1) > graythresh(sAMP);
    %figure;imshow(T(:,:,1) == 3);
    M = (T(:,:,1) == 3 | T(:,:,1) == 5) & ampMASK;
    %{
    skel = bwmorph(M,'thin',inf);
    fidx = find(M);
    newSIG = mean(g3(fidx,1:3));
    DIS = bsxfun(@minus,g3(:,1:3),newSIG);
    DIS = sum(DIS.*DIS,2).^.5;
    DIS = reshape(DIS,size(d1));
    RGB = cat(3,bindVec(T(:,:,1)),bindVec(T(:,:,2)),bindVec(T(:,:,3)));
    %}
    %{
    fidx = find(M);
    SAMP = zeros([numel(fidx) NF]);
    for slice = 1:size(T,3)
        tmp = AMP(:,:,slice);
        SAMP(:,slice) = tmp(fidx);
    end
    TARGET = mean(SAMP,1);
    TARGET = repmat(shiftdim(TARGET,-1),[size(AMP,1) size(AMP,2) 1]);
    DIST = sum((TARGET - AMP).^2,3).^.5;
    %}
    
    M = imresize(M,size(tmpR));
    %M = imerode(M,strel('disk',5));
    Re = regionprops(M,'MajorAxisLength','PixelIdxList');
    [J midx] = max([Re.MajorAxisLength]);
    M = zeros(size(M));
    M(Re(midx).PixelIdxList) = 1;
    M = imfill(M,'holes');
    M = imdilate(M,strel('disk',5));

    M = padarray(M, [R(2) R(2)]);
    M = imresize(M,size(tmpO));

    M = imclose(logical(M),strel('disk',50));
    skel = bwmorph(M,'thin',inf);
    skeleton = imdilate(skel,strel('disk',3));
    %{
    figure;
    imshow(M,[]);
    out = flattenMaskOverlay(tmpO,logical(skel));
    figure;
    imshow(out);
    %}
end