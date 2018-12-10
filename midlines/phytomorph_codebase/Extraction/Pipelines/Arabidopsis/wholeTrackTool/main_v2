
FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/straight/';
FileList = {};
FileExt = {'tiff','TIF'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);



for s = 1:20%:numel(SET)
    % read the first image from the stack
    I = imread(SET{s}{1});
    % track every Xth pixel
    pixelSKIP = 7;
    % 
    v = min(floor((size(I)-1)/pixelSKIP));
    I = I(1:(pixelSKIP*v+1),1:(pixelSKIP*v+1));
    [n1 n2] = ndgrid(1:pixelSKIP:size(I,1),1:pixelSKIP:size(I,2));
    bsz = factor(size(n1,1));
    % 
    targetBK = 25;
    [~,bidx] = min(abs(bsz - targetBK));
    blockSize = bsz(bidx);
    % create block size of 14 at X pixel spacing
    blockSize = 14;
   
    IDXy = im2col(n1,[blockSize blockSize],'distinct');
    IDXx = im2col(n2,[blockSize blockSize],'distinct');


    for e = 1:size(IDXx,2)
        IDX{e} = [IDXx(:,e) IDXy(:,e)];
    end


    [g1 g2] = ndgrid(1:size(I,1),1:size(I,2));


    F = IDX;
    func = cFlow('wholeTrack');
    func.setMCRversion('v840');
    func.setMemory(4000);
    for f = 1:numel(SET{s})
        [pth,nm,ext] = fileparts(SET{s}{f});
        tmpFile{f} = [nm ext];
        func.addLocalD(SET{s}{f});
    end
    pointList = {};
    %tmpFile = SET{s};
    pointListL = {};
    pointListH = {};
    SL = {};
    SE = {};
    for c = 1:numel(F)
        [pointListLs{s}{c} pointListEs{s}{c} SL{s}{c} SE{s}{c}] = func(tmpFile,F{c},0,[]);
        %[pointListL{c} pointListE{c} SL{c} SE{c}] = wholeTrack(tmpFile,F{c},0,miniSKIP);
    end

    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};

    func.submitDag(auth,700,700);
end

%%


    s = 5;
    TILE_L = [];
    TILE_E = [];
    VEL_L = [];
    VEL_E = [];
    cnt = 1;
    close all;
    for l = 1:numel(pointListLs{s})
        SZ{l} = [14 14];
        [UNITL] = loadComputeUnit(pointListLs{s}{l},SZ{l});
        [UNITE] = loadComputeUnit(pointListEs{s}{l},SZ{l});
        %{
        imshow(I,[]);
        hold on
        K1 = squeeze(UNITL(:,:,1,1));
        K2 = squeeze(UNITL(:,:,2,1));
        plot(K1,K2,'r.')
      
        J1 = squeeze(UNITL(:,:,1,end));
        J2 = squeeze(UNITL(:,:,2,end));
        plot(J1,J2,'b.')
        hold off
        drawnow
        %}
       
        TILE_L = [TILE_L,UNITL];
        TILE_E = [TILE_E,UNITE];
        if size(TILE_L,2) == ((size(n2,2)))
            %break
            cnt = 1;
            VEL_L = [VEL_L;TILE_L];
            VEL_E = [VEL_E;TILE_E];
            TILE_L = [];
            TILE_E = [];
        else
            cnt = cnt + 1;
        end
       
    end
    

    
    MAP_L = zeros(size(I,1),size(I,2),2,numel(SET{s}));
    for t = 1:size(VEL_L,4)
        MAP_L(:,:,1,t) = imresize(VEL_L(:,:,1,t),size(I),'bicubic');
        MAP_L(:,:,2,t) = imresize(VEL_L(:,:,2,t),size(I),'bicubic');
    end


    
    MAP_E = zeros(size(I,1),size(I,2),2,numel(SET{s}));
    for t = 1:size(VEL_E,4)
        MAP_E(:,:,1,t) = imresize(VEL_E(:,:,1,t),size(I),'bicubic');
        MAP_E(:,:,2,t) = imresize(VEL_E(:,:,2,t),size(I),'bicubic');
    end
    [g1 g2] = ndgrid(1:size(I,1),1:size(I,2));
    

    % mean of L flow
    FF = mean(sum(diff(MAP_L,1,4).^2,3).^.5,4);
    FF(FF > 5) = 0;
    imshow(FF,[]);

    [VEC] = makeItWork_ver0(MAP_L,I,.15,21);


    % dislay the final and first L
    SKIPP = 10;
    K1 = squeeze(VEL_L(1:SKIPP:end,1:SKIPP:end,1,1));
    K2 = squeeze(VEL_L(1:SKIPP:end,1:SKIPP:end,2,1));
    J1 = squeeze(VEL_L(1:SKIPP:end,1:SKIPP:end,1,end));
    J2 = squeeze(VEL_L(1:SKIPP:end,1:SKIPP:end,2,end));
    imshow(I,[]);hold on
    plot(K1(:),K2(:),'r.')
    plot(J1(:),J2(:),'b.')

    % display E movie
    SPEED = squeeze(sum(MAP_E.^2,3).^.5);
    for s = 1:size(SPEED,3)
        tmp = SPEED(:,:,s);
        tmp(tmp > 5) = 0;
        imshow(tmp,[0 5]);
        drawnow
    end

    
    STRAIN = [];
    for e = 1:size(SPEED,3)
        tmp = imfilter(SPEED(:,:,e),fspecial('gaussian',[71 71],31),'replicate');
        %imshow(tmp,[])
        %drawnow
        tmpS = [];
        for d = 1:size(MAP_E,3)
            tmpS(:,:,d) = imfilter(MAP_E(:,:,d,e),fspecial('gaussian',[21 21],7),'replicate');
        end
        [STRAIN(:,:,1,e) STRAIN(:,:,2,e)] = gradient(tmp);
        STRAIN(:,:,:,e) = STRAIN(:,:,:,e).*tmpS;
    end
    STRAIN = squeeze(sum(STRAIN,3));

    for e = 1:size(STRAIN,3)
        tmp = STRAIN(:,:,e);
        tmp = imfilter(tmp,fspecial('disk',21),'replicate');
        tmp(abs(tmp) > .03) = 0;
        imshow(tmp,[])
        drawnow
    end

    I = imread(SET{s}{1});
    pixelSKIP = 7;
    v = min(floor((size(I)-1)/pixelSKIP));
    I = I(1:(pixelSKIP*v+1),1:(pixelSKIP*v+1));

    figure, imshow(I);
    h = impoly;
    position = wait(h);
   
    MSK = poly2mask(position(:,1),position(:,2),size(I,2),size(I,1));
    [q1 q2] = impixel(I);

    [Z Z1] = sampleRegion(MAP_E,MAP_L,MSK,[q1 q2]);

    

    for t = 1:size(MAP,4)
        for d = 1:size(MAP,3)
            tmpZ = MAP(:,:,d,t);
            tmpZ(zidx) = 0;
            MAP(:,:,d,t) = tmpZ;
        end
    end

    dMAP = bsxfun(@minus,MAP,MAP(:,:,:,1));
    d1 = squeeze(dMAP(:,:,1,:));
    sz1 = size(d1);
    d1 = reshape(d1,[size(d1,1)*size(d1,2) size(d1,3)]);
    d2 = squeeze(dMAP(:,:,2,:));
    d2 = reshape(d2,[size(d2,1)*size(d2,2) size(d2,3)]);
    [S1 C1 U1 E1 L1 ERR1 LAM1] = PCA_FIT_FULL(d1,3);
    [S2 C2 U2 E2 L2 ERR2 LAM2] = PCA_FIT_FULL(d2,3);
    S1 = reshape(S1,[sz1(1:2) 1 sz1(3)]);
    S2 = reshape(S2,[sz1(1:2) 1 sz1(3)]);
    sMAP = cat(3,S1,S2);
    sMAP = bsxfun(@plus,sMAP,MAP(:,:,:,1));
    tmp = [];
    for d = 1
        tmp(:,:,d) = bindVec(reshape((C1(:,d)).*(C2(:,d)),[sz1(1:2) 1]));
    end
    imshow(tmp,[]);
    imshow(tmp < .1);

    imshow(I,[]);hold on;
    sv = fliplr(21:2:131);
    for sm = 1:numel(sv)
        s = imfilter(dM,fspecial('disk',sv(sm)));
        [~,idx] = max(s(:));
        [i1 i2] = ind2sub(size(s),idx);
        plot(i2,i1,'r*')
        drawnow
    end
    dMAP = [];
    for d = 1:size(sMAP,3)
        [dMAP(:,:,:,d,1) dMAP(:,:,:,d,2) dMAP(:,:,:,d,3)] = gradient(squeeze(sMAP(:,:,d,:)));
    end

    MSK = tmp < .75;

    [pt(1,1) pt(1,2) V] = impixel(I);
    F = imread(SET{1}{1});
    followFlow(F,MAP_L,pt)
    SKIPP = 5;
    DEL = VEL(:,:,:,end) - VEL(:,:,:,1);
    
    K1 = squeeze(VEL(1:SKIPP:end,1:SKIPP:end,1,1));
    K2 = squeeze(VEL(1:SKIPP:end,1:SKIPP:end,2,1));
    J1 = squeeze(VEL(1:SKIPP:end,1:SKIPP:end,1,end));
    J2 = squeeze(VEL(1:SKIPP:end,1:SKIPP:end,2,end));
    imshow(I,[]);hold on
    plot(K1(:),K2(:),'r.')
    plot(J1(:),J2(:),'b.')
    plot(MAP(pt(1,2),pt(1,1),1,1),MAP(pt(1,2),pt(1,1),2,1),'g*')
    plot(MAP(pt(1,2),pt(1,1),1,end),MAP(pt(1,2),pt(1,1),2,end),'k*')
 
    

end






%}