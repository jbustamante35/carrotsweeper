function [M] = hmm_ver0(ST)

    uST = mean(ST,4);
    [boundary currentArea MASK] = getKernelArea(uST,5000);
    [boundary currentArea iMASK] = getKernelArea(ST(:,:,:,1),5000);
    [boundary currentArea fMASK] = getKernelArea(ST(:,:,:,end),5000);
    E = edge(MASK);
    dB = bwdist(E);
    eMASK = dB < 40;
    eMASK = ((eMASK - iMASK) == 1);
    %eMASK = imdilate(eMASK,strel('disk',3));
    eMASK = (fMASK - iMASK)==1;
    fidx = find(eMASK);
    
    sz = size(ST);
    SIG = reshape(ST,[prod(sz(1:2)) sz(3) sz(4)]);
    SIG = permute(SIG,[2 3 1]);
    SIG = SIG(:,:,fidx);
    
    FGK = [];
    BGK = [];
    SKIP = 1;
    iVec =  1:SKIP:50;
    for e = 1:numel(iVec)
        tmpI = ST(:,:,:,iVec(e));
        [boundary currentArea MASK] = getKernelArea(tmpI,5000);
        MASK = imerode(MASK,strel('disk',5));
        fidx1 = find(MASK==1);
        fidx0 = find(MASK==0);
        tmpI = reshape(tmpI,[size(tmpI,1)*size(tmpI,2) size(tmpI,3)]);
        FGK = [FGK;tmpI(fidx1(1:1:end),:)];
        BGK = [BGK;tmpI(fidx0(1:1:end),:)];
        e
        numel(iVec)
    end
    FGK = double(FGK);
    BGK = double(BGK);
    background_node = hmm_node('background');
    foreground_node = hmm_node('foreground');
    bk_to_bk = constantTransitionFunction(.90);
    bk_to_fg = constantTransitionFunction(.10);
    fg_to_fg = constantTransitionFunction(1);
    fg_to_bk = constantTransitionFunction(0);
    background_node.attachNode(background_node,bk_to_bk);
    foreground_node.attachNode(foreground_node,fg_to_fg);
    background_node.attachNode(foreground_node,bk_to_fg);

    fg_distribution = myProb(mean(FGK,1),cov(FGK));
    bg_distribution = myProb(mean(BGK,1),cov(BGK));
    fg_distribution.fitToGMM(FGK,2,100);
    bg_distribution.fitToGMM(BGK,2,100);
    foreground_node.attachDistribution(fg_distribution,1);
    background_node.attachDistribution(bg_distribution,1);
    hmm = my_hmm();
    hmm.addNode(background_node);
    hmm.addNode(foreground_node);
    hmm.dn = 1;
    ST = double(ST);
    %{
    [J BOXP] = imcrop(ST(:,:,:,1)/255);
    
    
    rST = [];
    cnt = 1;
    SKIP = 5;
    for e = 1:SKIP:size(ST,4)
        %rST(:,:,:,e) = imresize(ST(:,:,:,e),.05);
        rST(:,:,:,cnt) = imcrop(ST(:,:,:,e),BOXP);
        cnt = cnt + 1;
    end
    rST = double(rST);
    
    SIG = permute(rST,[1 2 3 4]);
    sz = size(SIG);
    SIG = reshape(SIG,[prod(sz(1:2)) sz(3) sz(4)]);
    SIG = permute(SIG,[2 3 1]);
    %}
    seg = [];
    SEG = {};
    SIG = double(SIG);
    
    % PRE
    % compute the log determinant of covariance
    
    
    seg = [];
    e=1;
    tmp = SIG(:,:,e);
    LAB = ones(3,1);
    %LAB = 1;
    junk = hmm.Viterbi(tmp,LAB);
    parfor e = 1:size(SIG,3)
        tic
        tmp = SIG(:,:,e);
        LAB = ones(3,1);
        %LAB = 1;
        seg(e,:) = hmm.Viterbi(tmp,LAB);
        SEG{e} = [seg(e,:);tmp];
        toc*size(SIG,3)/60/12
        e
        size(SIG,3)
    end
    
    
    hmm.update(SEG,ones(3,1));
    hmm.needPRE = 1;
    
    % make new masks    
    nMASK = [];
    for e = 1:size(ST,4)
        tmp = iMASK;
        tmp(fidx) = (seg(:,e) == 2);
        %tmp = bwareaopen(tmp,5000);
        nMASK(:,:,e) = tmp;
        
    end
    
    close all
    AREA = [];
    for tm = 1:size(nMASK,3)
        tmp = double(nMASK(:,:,tm));
        [boundary currentArea MASK2] = getKernelArea(ST(:,:,:,tm),5000);
        tmp = imclose(tmp,strel('disk',3));
        tmp = imfilter(tmp,fspecial('average',[11 11]));
        tmp = tmp > graythresh(tmp);
        AREA(tm) = sum(tmp(:));
        AREA2(tm ) =sum(MASK2(:));
        %tmp = iMASK;
        B = bwboundaries(tmp);
        DB{tm} = B{1};
        imshow(ST(:,:,:,tm)/255,[]);
        hold on
        plot(B{1}(:,2),B{1}(:,1),'r')
        if ~isempty(boundary)
            plot(boundary(:,2),boundary(:,1),'b')
        end
        hold off
        title(num2str(tm))
        drawnow
    end
    
    
    for e = 1:numel(DB)
        d = size(DB{e},1);
       plot3(DB{e}(:,1),DB{e}(:,2),e*ones(d,1));
       hold on
    end
    % view seg
    seg = reshape(seg,[size(rST,1) size(rST,2) size(rST,4)]);

    
    for tm = 1:size(seg,3)
        imshow(seg(:,:,tm)==2,[]);
        drawnow
    end
    
    
    
    % try the expanding boundary
    [boundary currentArea iMASK] = getKernelArea(ST(:,:,:,1),5000);
    
     % make new masks
    for e = 1:size(ST,4)
        tmp = iMASK;
        nMASK(:,:,e) = tmp;
    end
    sz = size(ST);
    SIG = reshape(ST,[prod(sz(1:2)) sz(3) sz(4)]);
    SIG = permute(SIG,[2 3 1]);
    eMASK = edge(iMASK);
    fidx = find(eMASK);
    flag = 1;
    while flag
        
        tSIG = SIG(:,:,fidx);
        seg = [];
        parfor e = 1:size(tSIG,3)
            tic
            tmp = tSIG(:,:,e);
            LAB = ones(3,1);
            %LAB = 1;
            seg(e,:) = hmm.Viterbi(tmp,LAB);
            toc*size(tSIG,3)/60/10
            e
            size(tSIG,3)
        end
        % make new masks
        for e = 1:size(ST,4)
            tmp = nMASK(:,:,e);
            tmp(fidx) = (seg(:,e) == 2);
            %tmp = bwareaopen(tmp,5000);
            nMASK(:,:,e) = tmp;
        end
        flag = ~all(seg(:,end)==1);
        sum(seg(:,end)==1)
        iMASK = imdilate(iMASK,strel('disk',1));
        eMASK = edge(iMASK);
        fidx = find(eMASK);
    end
end
