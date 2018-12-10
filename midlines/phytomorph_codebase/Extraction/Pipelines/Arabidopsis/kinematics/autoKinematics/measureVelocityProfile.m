function [domain,vP] = measureVelocityProfile(pointList,DST,DSS,np)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % step 2: center the tip
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xm = pointList;        
    for t = 1:size(Xm,3)
        Xm(:,:,t) = Xm(:,:,t) - repmat(Xm(1,:,t),[size(Xm(:,:,t),1) 1 1]);
    end
    tmpX = Xm(1:DSS:end,:,1:DST:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the position values
    [J,POS,VEL] = gradient(tmpX);
    % POS = diff(tmpX,1,1); 
    % POS = gradient(tmpX);
    POS = sum(POS.*POS,2).^.5;
    h1 = mean(POS(:));
    %POS = cat(1,zeros(1,1,size(POS,3)),POS);
    POS = cumsum(POS,1);
    POS = squeeze(POS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the velocity values
    VEL = sum(VEL.*VEL,2).^.5;
    h2 = mean(VEL(:));
    VEL = squeeze(VEL);
    VEL = VEL * DSS^-1;
    TM_DOMAIN = repmat(1:size(POS,2),[size(POS,1) 1]);
    NP = (max(POS(:,1)) - min(POS(:,1)))/np;
    [iD1 iD2] = ndgrid(linspace(min(POS(:)),max(POS(:,1)),round(NP)),1:size(POS,2));
    h1 = diff(iD1,1,1);
    h1 = mean(h1(:));
    h2 = diff(iD2,1,2);
    h2 = mean(h2(:));
    
    vP = griddata(POS/h1,TM_DOMAIN/h2,VEL,iD1*h1.^-1,iD2*h2.^-1);
    domain = cat(3,iD1,iD2);
end