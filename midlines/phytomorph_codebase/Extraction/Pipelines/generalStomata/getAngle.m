function [result] = getAngle(I,P,net,Z0,trainedNetwork_Skew_local,ZS,netMAJOR,Z1,netMINOR,Z2,chop,grader)
    
    % generate n-sub index
    [xPos(2),xPos(1)] = ind2sub(size(I),P);
    patch = I((xPos(2)-chop):(xPos(2)+chop),(xPos(1)-chop):(xPos(1)+chop),:);
    % displace the domain
    anglePatch = zscore(patch(:));
    anglePatch = reshape(anglePatch,size(patch));
    ang = net.predict(anglePatch).*Z0(2,:) + Z0(1,:);
    ang = ang(end);
    % sample the square grid
    [x1,x2] = ndgrid(linspace(-chop,chop,2*chop+1),linspace(-chop,chop,2*chop+1));
    xD = [x2(:) x1(:)];
    PARA(1) = ang;
    PARA(2:3) = 1;
    PARA(4:5) = 0;
    

    %{
        figure;
        imshow(patch,[]);
        hold on
        plot((size(patch,1)-1)/2+1,(size(patch,2)-1)/2+1,'.')

        [T] = generateDomainTransformation(PARA);
        T(:,3) = T(:,3) + [31;31];
        TH = linspace(-pi,pi,100);
        dX10 = [cos(TH)' sin(TH)'];
        [dX1] = transformDomain(dX10,T);
        plot(dX1(:,1),dX1(:,2),'g');


        imshow(patch,[]);
        cP = reshape(yP(1,1:8),[2 4])';
        hold on
        plot(cP(:,1)+31,cP(:,2)+31,'r.')
        plot(31,31,'g*');
        gM = norm(cP(1,:)-cP(3,:))/2;
        gu = norm(cP(2,:)-cP(4,:))/2;
        Mvec = -(cP(3,:)-cP(1,:));
        uvec = -(cP(4,:)-cP(2,:));
        ang(1) = yROT(eP,end);
        ang(2) = atan2(Mvec(2),Mvec(1));
        ang(3) = atan2(uvec(2),uvec(1)) - pi/2;
        paraDisplay = [mean(ang) gM gu 0 0];
        [Tx] = generateDomainTransformation([paraDisplay 0 0]);
        TH = linspace(-pi,pi,50);
        dispX = [cos(TH)' sin(TH)' ones(50,1)];
        dispX = (Tx*dispX')';
        plot(dispX(:,1)+31,dispX(:,2)+31,'b')
        plot(cP(1,1)+31,cP(1,2)+31,'m*')



    %}
    
    [T] = generateDomainTransformation(PARA);
    [xD] = transformDomain(xD,T);
    xD = double(xD);
    [xI] = myInterp2Sampler(I,P,xD,size(x1));
    


    fsz = size(xI);
    xI = zscore(xI(:));
    xI = reshape(xI,fsz);
    
    majorAxis = netMAJOR.predict(xI).*Z1(2) + Z1(1);
    %minorAxis = netMINOR.predict(xI).*Z2(2) + Z2(1);
    PARA(2) = majorAxis(1);
    PARA(3) = majorAxis(2);


    % sample the disk
    [theta rho] = ndgrid(linspace(-pi,pi,100),linspace(0,1.5,50));
    xD = [rho(:).*cos(theta(:)) rho(:).*sin(theta(:))];

    [x1,x2] = ndgrid(linspace(-1.5,1.5,50),linspace(-1.5,1.5,50));
    %[x1,x2] = ndgrid(linspace(-3,3,100),linspace(-3,3,100));
    xD = [x2(:),x1(:)];
    
    [T] = generateDomainTransformation(PARA);
    %T(:,3) = T(:,3) + xPos';
    [xD] = transformDomain(xD,T);
    xD = double(xD);
    [xI] = myInterp2Sampler(I,P,xD,size(x1));


    fsz = size(xI);
    fxI = zscore(xI(:));
    xI = reshape(xI,fsz);

    fxI = xI;
    %{
    for loop = 1:3
        % find the correction angle
        skewAngle = trainedNetwork_Skew_local.predict(xI)*ZS(2) + ZS(1);
        x = cos(skewAngle);
        y = sin(skewAngle);
        y = (PARA(3).*PARA(2).^-1).*y;
        correctionAngle = atan2(y,x);
        PARA(1) = PARA(1) + correctionAngle;


        % resample
        xD = [x2(:),x1(:)];
        [T] = generateDomainTransformation(PARA);
        %T(:,3) = T(:,3) + xPos';
        [xD] = transformDomain(xD,T);
        xD = double(xD);
        [xI] = myInterp2Sampler(I,P,xD,size(x1));


        fsz = size(xI);
        fxI = zscore(xI(:));
        xI = reshape(xI,fsz);

    end


    % re-check major and minor
    % sample the square grid
    [x1,x2] = ndgrid(linspace(-chop,chop,2*chop+1),linspace(-chop,chop,2*chop+1));
    xD = [x2(:) x1(:)];
    PARA(2:3) = 1;
    [T] = generateDomainTransformation(PARA);
    [xD] = transformDomain(xD,T);
    xD = double(xD);
    [fxI] = myInterp2Sampler(I,P,xD,size(x1));
    majorAxis = netMAJOR.predict(fxI).*Z1(2) + Z1(1);
    minorAxis = netMINOR.predict(fxI).*Z2(2) + Z2(1);
    PARA(2) = majorAxis(1);
    PARA(3) = minorAxis(2);


    fsz = size(fxI);
    fxI = zscore(fxI(:));
    fxI = reshape(fxI,fsz);


    %}


    xI = fxI(:);
    if ~isempty(grader)
           grade = grader.predict(reshape(xI,[50 50 1]));
    else
        grade = NaN;
    end

    

    
    PARA(4:5) = [];
    PARA = PARA';
    result = [xI' PARA(:)' grade];
    %PARA = xI;
    
end