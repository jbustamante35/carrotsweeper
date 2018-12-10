function [T,rawI] = ixExtractor(rawI,toOp,trimValues,pointSet,opDomain,domainSize,opFunc,disp)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load is need for raw image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(rawI)
        rawI = double(imread(rawI));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load is need for raw image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(toOp)
        % if nms file else tiff - left over from pre nms-reader
        toOp = double(imread(toOp));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    mag = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resize the image and strip edges
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mag ~= 1
        rawI = imresize(rawI,mag);
        toOp = imresize(toOp,mag);
    end
    for r = 1:4
        rawI(:,1:trimValues(2)) = [];
        rawI = imrotate(rawI,90);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create sample disk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
    X = n1.*cos(n2);
    Y = n1.*sin(n2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create sample point(s) grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(P)
        [d1,d2] = ndgrid((R(2)+1):(size(toOp,1)-R(2)),(R(2)+1):(size(toOp,2)-R(2)));
        PZ = size(d1);
        P = [d1(:),d2(:)];
    elseif numel(P) == 1
        [P(:,2),P(:,1),~] = impixel(toOp);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    
    
    
    
    %{
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init store
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% init retVec1
    numFeatures1 = 2;
    pointsVec1 = PZ;
    numPerFeature1 = [N(1) NF];
    numSlice1 = size(toOp,3);
    retVec1 = zeros(prod(pointsVec1),prod(numPerFeature1)*numFeatures1*numSlice1);
    %%%% init retVec2
    numFeatures2 = 1;
    pointsVec2 = PZ;
    
    %numPerFeature2 = [N(1) NF NF];
    % 
    numPerFeature2 = [N(1) ((NF^2)-NF)/2];
    %numPerFeature2 = [N(1) NF];
    % time based corr
    numPerFeature2 = [N(1)-1 NF];
    numSlice2 = size(toOp,3);
    retVec2 = zeros(prod(pointsVec2),prod(numPerFeature2)*numFeatures2*numSlice2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % peform extraction and analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['**********************************\n']);
    fprintf(['Starting opFunc over pointSet:']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each slice
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for slice = 1:size(toOp,3)
        
        %tmpOp = toOp(:,:,slice);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each point
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['start:Whole extracting bug eye image(s)\n']);
        str = clock;
        [subI] = bugEyeViaT(toOp,pointSet(1,:),opDomain,domainSize);
        subI = zeros(size(subI,1),size(pointSet,1));
        parfor p = 1:size(pointSet,1)
            fprintf(['start:Extracting bug eye image ' num2str(p) ':' num2str(size(pointSet,1)) '.\n']);tic
            subI(:,p) = bugEyeViaT(toOp,pointSet(p,:),opDomain,domainSize);
            eTime = toc;
            fprintf(['end:Extracting bug eye image ' num2str(p) ':' num2str(size(pointSet,1))  '@' num2str(eTime) '.\n']);
            rTime = mean(eTime)*(size(pointSet,1)-p);
            fprintf(['esti: time remaining : ' num2str(rTime) '\n']);
        end
        fprintf(['end:Whole extracting bug eye image(s) ' num2str(etime(clock,str)) '\n']);
        
        
            
        T = opFunc(subI);
        %T = newF;
        %T = glueFrozenTensors(subI,newF);
        
        
        
            
            %[rZ,rDIS,AMP,ANG] = robustGetOmega(P(p,:),toOp,R,N,NF,n,dK,ikm,false);
            
            
            
            
            
            %rDIS = permute(rDIS,[3 1 2]);
            %{
            tic
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get location
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %LOC = [d1(p),d2(p)];
            LOC = P(p,:);
            Xp = X + LOC(2);
            Yp = Y + LOC(1);

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % display the sample scan
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if disp
                imshow(I);
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % perform interpolation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            F = ba_interp2(tmpOp,Xp,Yp,'cubic');
            %F = bsxfun(@minus,F,mean(F,2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % perform analysis
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fprintf(['Starting fourier calculation.\n']);
            fF = fft(F,[],2);
            %fprintf(['Ending fourier calculation.\n']);
            
            
            %fprintf(['Starting angle extraction.\n']);
            ang = angle(fF(1:N(1),1:NF));  
          
            %if nargout == 4
               % ORG(:,:,p) = F;
            %end
            
            
            [DIS] = calcDM(ang);
            toc
            %}
            %ang = unwrap(ang,[],1);
            %fprintf(['Ending angle extraction.\n']);
            
            %{
            %fprintf(['Starting wave disperion calculation.\n']);
            targetFrequency = 2*pi.*(0:(NF-1)).^-1;
            minJump = [];
            for rho = 1:size(ang,1)
                for dispersion = 1:size(ang,2)
                    baseAngle = ang(rho,dispersion);
                    [~,minJump(dispersion,:,rho)] = phD(baseAngle,ang,targetFrequency);
                end
            end
            minJump(isnan(minJump(:))) = 0;
            minJump = permute(minJump,[3 1 2]);
            %fprintf(['Ending wave disperion calculation.\n']);
            %}
            
            %fprintf(['Starting amplitude extraction.\n']);
            %fF = abs(fF(1:N(1),1:NF));
            %fprintf(['Ending amplitude extraction.\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % store the data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %retVec1(p,:) = [fF(:);ang(:)]';
            %retVec1(p,:) = [AMP(:);ANG(:)]';
            %retVec1(p,:) = fF(:)';
            %retVec2(p,:) = rDIS(:)';
            %retVec2(p,:) = ang(:)';


            
            
        %end
    %end
    fprintf('\n');
    fprintf(['Ending FFT:\n']);
    fprintf(['**********************************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reshape output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %retVec1 = reshape(retVec1,[pointsVec1,numPerFeature1,numFeatures1,numSlice1]);
    %retVec2 = reshape(retVec2,[pointsVec2,numPerFeature2,numFeatures2,numSlice2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
