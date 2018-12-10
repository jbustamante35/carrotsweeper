function [BUG,tmpR] = buyEye_NET(I,R,N,mag,disp,funcIN,funcOUT)    
    % I: = image to sample
    % R: = radius vector of inner to outer disk
    % N: = number of points
    % mag:= resize the image
    % func:= to apply to the sample
    
    if isdeployed
       % parpool(1)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resize the image
    tmpI = imresize(I,mag);
    
    if ~isempty(funcIN)
        tmpI = funcIN(tmpI);
    end
    
    tmpR = tmpI;
    for r = 1:4
        tmpR(:,1:R(2),:) = [];
        tmpR = imrotate(tmpR,90);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create sample disk
    [n1 n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
    [d1 d2] = ndgrid(R(2):(size(tmpI,1)-R(2)+1),R(2):(size(tmpI,2)-R(2)+1));
    X = n1.*cos(n2);
    Y = n1.*sin(n2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['**********************************\n']);
    fprintf(['Starting BUG NET:']);
    
    if ~isempty(funcOUT)
        %% try one point passed through the function
        p=1;
        LOC = [d1(p) d2(p)];
        Xp = X + LOC(2);
        Yp = Y + LOC(1);
        % display the sample scan
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
        FF = applyScan(F,funcOUT);
        sz3 = size(FF);
    else
        sz3 = [size(X) size(tmpR,3)];
    end
    
    BUG = zeros([numel(d1) sz3]);
    % for each point
    parfor p = 1:numel(d1)
        LOC = [d1(p) d2(p)];
        Xp = X + LOC(2);
        Yp = Y + LOC(1);
        % display the sample scan
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
        
        if ~isempty(funcOUT)
            F = applyScan(F,funcOUT);
        end
        
        BUG(p,:,:,:) = F;
       
    end
    BUG = reshape(BUG,[size(d1) sz3]);
    
    fprintf('\n');
    fprintf(['Ending bug NET:\n']);
    fprintf(['**********************************']);
end