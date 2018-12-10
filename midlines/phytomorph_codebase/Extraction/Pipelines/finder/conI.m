function [vec] = conI(I,reSize,gridSize,func)
    for r = 1:numel(reSize)
        rI = imresize(I,reSize(r));
        para.patchSize = gridSize;
        para.func = func;
        
        vec(r) = func1(rI,para);
    end
    vec = sum(vec);
    
end

%{
   

    % get other image
    I = double(imread('/home/nate/Downloads/000000.TIF'));



    % get video image
    videoR = VideoReader('/home/nate/Downloads/GOPR0495.MP4');
    nFrames = videoR.NumberOfFrames;
    I = double(readFrame(videoR));

    % read in frames
    [mov] = readMP4('/home/nate/Downloads/GOPR0495.MP4',3);
    [mov2] = readMP4('/home/nate/Downloads/GOPR0496.MP4',3);
    [mov3] = readMP4('/home/nate/Downloads/GOPR0500.MP4',3);

    % analysis of std of movie
    s = [];
    for e = 1:size(mov3,4)
        g = rgb2gray(mov3(:,:,:,e));
        g = stdfilt(g);
        s(e) = mean(g(:));
        e
    end


    % first movie stack
    mv = mov(:,:,:,[40:200 234:size(mov,4)]);
    target = [ones(200-40+1,1);zeros(size(mov,4)-234+1,1)];
    mv = double(mv)/255;

    % second movie stack
    mv2 = mov2(:,:,:,[40:120 164:size(mov2,4)]);
    target2 = [ones(120-40+1,1);zeros(size(mov2,4)-164+1,1)];
    mv2 = double(mv2)/255;


    % second movie stack
    mv3 = mov3(:,:,:,[50:135 150:size(mov3,4)]);
    target3 = [ones(135-50+1,1);zeros(size(mov3,4)-150+1,1)];
    mv3 = double(mv3)/255;

    blockSize = [25 25];
    [vec] = prep2(mv,blockSize,linspace(.1,.8,8));
    [vec2] = prep2(mv2,blockSize,linspace(.1,.8,8));
    [vec3] = prep2(mv3,blockSize,linspace(.1,.8,8));
    X = [vec;vec2;vec3];
    Y = [target;target2;target3];

    for m = 1:50
        for h = 1:size(X,1)
            idx = setdiff(1:size(X,1),h);
            [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(X(idx,:),Y(idx),m);
            sc = [ones(1,1) X(h,:)]*BETA;
            tmp(h) = sc'*Y(h);
        end
        score(m) = mean(tmp);
        m
        plot(score)
        drawnow
    end


    for m = 1:50
        [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(X,Y,m);
        lam = reshape(BETA(2:end),[blockSize 3]);
        imshow(imresize(cat(3,bindVec(lam(:,:,1)),bindVec(lam(:,:,2)),bindVec(lam(:,:,3))),4),[])
        drawnow
    end

    [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(X,Y,10);
    sc = [ones(size(X,1),1) X]*BETA;




    
    matchScale(mv(:,:,:,30),BETA,blockSize,linspace(.1,.85,9));
    matchScale(mv2(:,:,:,50),BETA,blockSize,linspace(.1,1,9),E,U);
    matchScale(mv3(:,:,:,50),BETA,blockSize,linspace(.1,1,9));
    


    reSize = linspace(.2,.5,3);
    gridSize = [21 21 3];
    vec = rand(gridSize);
    vec = reshape(vec,[21*21 3]);
    for e = 1:size(vec,2)
        vec(:,e) = vec(:,e) / norm(vec(:,e));
    end
    
    
    
    func = @(m,k)myOp1(m,vec',k);      
    v = conI(I,reSize,gridSize(1:2),func);
    

    vec = rand(gridSize);
    vec = reshape(vec,[21*21 3]);
    for e = 1:size(vec,2)
        vec(:,e) = vec(:,e) / norm(vec(:,e));
    end
    reSize = linspace(.2,.3,3);
    gridSize = [21 21 3];
    tFunc = @(x)trainPuppy(x,mv(:,:,:,1:6:end),target(1:6:end));
    opt = optimoptions('fminunc','Display','iter');
    x = fminunc(tFunc,vec(:),opt);
    


    reSize = linspace(.2,.5,3);
    BLOCK = prep(mv(:,:,:,:),reSize,gridSize);
    sz = size(BLOCK);
    B = reshape(BLOCK,[prod(sz(1:2)) prod(sz(3))]);
    T = target;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(B',T,5);
   

    [vec] = prep2(mov,[21 21],linspace(.1,.5,5));
    tvec = vec([40:200 234:size(mov,4)],:);
    [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(tvec,target,7);
    sc = [ones(size(tvec,1),1) tvec]*BETA;


    matchScale(mv(:,:,:,30),BETA,[21 21],linspace(.1,.5,5));

    B = squeeze(sum(BLOCK,2));
    sz = size(B);
    B = reshape(B,[prod(sz(1:2)) prod(sz(3))]);
    T = target(1:3:end);
    [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(B',T,5);





%}