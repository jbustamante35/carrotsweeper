function [returnVector] = generatePARAfromNetworks_ver2(I,P,xD,sz,trainedNetworks,Z0,maxIter,gradeRatio,simplexPara,mag,disp,grader,imageGradFunc,statEnergyMix)




   
    [returnVector] = generatePARAfromNetworks_ver3(I,P,xD,sz,trainedNetworks,Z0,maxIter,false);

    
    returnGrade = NaN*ones(size(returnVector,1),1);


    %[M1,M2] = ndgrid(linspace(-1.5,1.5,50),linspace(-1.5,1.5,50));
    [M1,M2] = ndgrid(linspace(-2,2,100),linspace(-2,2,100));
    mD = [M2(:),M1(:)];
    
    %{
    [M1 M2] = ndgrid(linspace(0,2,50),linspace(-pi,pi,200));
    mD = [M1(:).*cos(M2(:)) M1(:).*sin(M2(:))];
    %}

    [H1,H2] = ndgrid(linspace(-30,30,61),linspace(-30,30,61));
    hD = [H2(:),H1(:)];
    %{
    [H1 H2] = ndgrid(linspace(0,30,31),linspace(-pi,pi,100));
    hD = [H1(:).*cos(H2(:)) H1(:).*sin(H2(:))];
    %}

    % view patch deformed
    [V1,V2] = ndgrid(linspace(-2,2,100),linspace(-2,2,100));
    vD = [V2(:),V1(:)];


    [P(1,2),P(1,1)] = ind2sub(size(I),P);
    %mag = [1 .5 .25];
   
    returnVector(1,4:5) = P;
    returnVector(2,4:5) = P; 


    geoReturnFunc = @(Q)generatePARAfromNetworks_ver3(I,Q,xD,sz,trainedNetworks,Z0,1,false);

    %{
    % grade with local diff
    initPara = returnVector(2,:);
    startPara = initPara;
    startPara(4:5) = P(1,:);
    %func = @(X)stomataMinObjectiveFunc(I,mD,size(M1),grader,P(1,:),[1 2],zeros(size(initPara)),geoReturnFunc,X);
    func = @(X)stomataMinObjectiveFunc(I,mD,size(M1),grader,P(1,:),P,[1 2],zeros(size(initPara)),geoReturnFunc,gradeRatio,disp,startPara,X);
               
    options = optimoptions('fminunc','FiniteDifferenceType','central','TypicalX',[pi/2 20 10 5 5],'FiniteDifferenceStepSize',[.0028 1/20 1/10 2/5 2/5]);
    iterPara = fminunc(func,initPara,options);
    returnVector = [returnVector;iterPara];

    iterP = P(1,:) + iterPara(4:5);
    iterPara(4:5) = 0;
    %}
    % store iterPara



    % simplex search
    gradeVecS = {[1 2],[3 4],[5 6]};
    gradeVecS = {[1:6]};
    gradeVecS = {[1 2],[1],[5 6]};
    %gradeVecS = {[7:8]};
    %gradeVecS = {[1:2 7]};
    gradeVecS = {[1]};
    options = optimset('TolFun',.1,'TolX',.1);
    for gr = 1:numel(gradeVecS)

        initPara = returnVector(1,:);
        initPara(4:5) = P(1,:);

        for m = 1:numel(mag)
        %{
            baseLine = mag(m)*simplexPara;
            deltaLB = [pi/4 5 5 5 5];
            deltaUB = [pi/4 5 5 5 5];
            startPara = initPara;
            fixerPara = baseLine - initPara;
            LB = baseLine - deltaLB;
            UB = baseLine + deltaUB;
            LB = initPara - deltaLB;
            UB = initPara + deltaUB;
            LB(3) = max(LB(3),5);
        
            fixerPara = 0;
            func = @(X)stomataMinObjectiveFunc(I,mD,size(M1),grader,[0 0],P,gradeVecS{gr},fixerPara,geoReturnFunc,gradeRatio,disp,startPara,X);
            options = optimoptions('particleswarm','SwarmSize',1000,'useParallel',true,'Display','iter','SocialAdjustmentWeight',10);
            x = particleswarm(func,5,LB,UB,options);
%}

            %baseLine = mag(m)*[20*pi/4 200 100 100 100];
            baseLine = mag(m)*simplexPara;
            startPara = initPara;
            fixerPara = baseLine - initPara;
            func = @(X)stomataMinObjectiveFunc(I,vD,size(V1),hD,size(H1),mD,size(M1),grader,[0 0],P,gradeVecS{gr},fixerPara,geoReturnFunc,gradeRatio,disp,startPara,X,imageGradFunc,statEnergyMix);
            [iterParaN,fval,exitflag,output] = fminsearch(func,baseLine,options);
            
            [minE,particleList] = func(iterParaN);

            fval = repmat(fval,[size(particleList,1) 1]);


            iterParaN = iterParaN - fixerPara;
            initPara = iterParaN;


            iterPN = iterParaN(4:5);
            iterParaN(4:5) = 0;
            
            toRet = particleList;
            

        end
        disP(gr).P = iterPN;
        disP(gr).para = iterParaN;

        
        returnVector = [returnVector;toRet];
        returnGrade = [returnGrade;fval];
    end
    

    displayPoint = P(1,:);
    displayPara = returnVector(2,:);

    %{
    close all

    dP = x(4:5);
    x(4:5) = 0;
    displayStoma(I,x,dP,'','r',1);
    displayStoma(I,displayPara,displayPoint,'','r',1);
    % display local diff method
    %displayStoma('',iterPara,iterP,'','g',1);
    CL = {'b' 'c' 'y'};
    for gr = 1:numel(disP)
        displayStoma('',disP(gr).para,disP(gr).P,'',CL{gr},1);
    end
   
    %}
    returnVector = [returnVector returnGrade];

    %waitforbuttonpress



end