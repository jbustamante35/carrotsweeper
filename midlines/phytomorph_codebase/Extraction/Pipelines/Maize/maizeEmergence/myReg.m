function [] = myReg(refI,movI,uCL,covCL)

        %refI = imresize(refI,.5);
        %movI = imresize(movI,.5);
        %movI = circshift(movI,[4 4]);

        sz = size(refI);
        refI_CL = reshape(refI,[prod(sz(1:2)) sz(3)]);
        refP = mvnpdf(refI_CL,uCL,covCL);
        refP = reshape(refP,sz(1:2));


        sz = size(movI);
        movI_CL = reshape(movI,[prod(sz(1:2)) sz(3)]);
        movP = mvnpdf(movI_CL,uCL,covCL);
        movP = reshape(movP,sz(1:2));


        movM = -log(movP) < -2;
        movM = bwareaopen(movM,500);
        movM = imfill(movM,'holes');
        movD = bwdist(movM) - bwdist(~movM);


        refM = -log(refP) < -2;
        refM = bwareaopen(refM,500);
        refM = imfill(refM,'holes');
        refD = bwdist(refM)  - bwdist(~refM);
            
        refM = refD < 50 & refD > 0; 
        movM = movD < 50 & movD > 0;
        

        masterM = refM & movM;


        fidx = find(~masterM);

        nref = refI;
        nmov = movI;
        for e = 1:3
            tmpI = nref(:,:,e);
            tmpI(fidx) = NaN;
            nref(:,:,e) = tmpI;


            tmpI = nmov(:,:,e);
            tmpI(fidx) = NaN;
            nmov(:,:,e) = tmpI;
        end


        R = regionprops(logical(bwlarge(imdilate(refM,strel('disk',5,0)))),'BoundingBox');
        refI = imcrop(refI,R(1).BoundingBox);
        movI = imcrop(movI,R(1).BoundingBox);
        func = @(X)myRegMeasure(rgb2gray(refI),refM,rgb2gray(movI),movM,X);
  
    nvars = 3 ;
    LB = [-5,-10,-10] ;
    UB = [5,10,10] ;
    options.PopInitRange = [[-5;5],[-5;5],[-5;5]] ;


 [xOpt,fval,exitflag,output,population,scores] =  pso(func,3,[],[],[],[],[],[],[],options,[0 0 0]);

        %T = imregcorr(rgb2gray(movI).*movM, rgb2gray(refI).*refM, 'rigid');
        %[output] = dftregistration(fft2(refM),fft2(movM),2);

        options = optimoptions(@fminunc,'Display','iter','DiffMinChange',.1,'FiniteDifferenceType','central','UseParallel',true);
        func = @(X)myRegMeasure(rgb2gray(refI),refM,rgb2gray(movI),movM,X);
        options = optimset('Display','iter','UseParallel',true);
        x = fminsearch(func,[0 0 0],options);
        x = fminunc(func,[0 0 0],options);

       

        T = [[cos(para(1)) sin(para(1)) 0];[-sin(para(1)) cos(para(1)) 0];[para(2) para(3) 1]]
        T = affine2d(T);

        






end