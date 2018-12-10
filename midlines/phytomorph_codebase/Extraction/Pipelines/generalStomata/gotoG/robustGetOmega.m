function [Z,DIS,AMP,ANG] = robustGetOmega(P,I,R,N,NF,n,dK,ikm,disp)
    delta = linspace(-pi,pi,n);
    

    [n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
    for r = 1:numel(delta)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create sample disk
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X(:,r) = n1(:).*cos(n2(:)-delta(r));
        Y(:,r) = n1(:).*sin(n2(:)-delta(r));
        %X(:,r) = n1(:).*cos(n2(:));
        %Y(:,r) = n1(:).*sin(n2(:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    

    p = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get location
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xp = X + P(p,2);
    Yp = Y + P(p,1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n > 1
        F = zeros(size(Xp,1),size(delta,2));
        for r = 1:numel(delta)
            %tmpI = imrotate(I,delta(r)*180/pi,'bicubic','crop');
            %GOUT{r} = tmpI;
            F(:,r) = ba_interp2(I,Xp(:,r),Yp(:,r),'cubic');
        end
    else
        F = ba_interp2(I,Xp,Yp,'cubic');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    F = reshape(F,[size(n1) n]);
    
    
    if disp
        for e = 1:size(F,3)
            imshow(F(:,:,e),[]);
            drawnow
        end
    end
    F = permute(F,[2 1 3]);
    fsz = size(F);
    if numel(fsz) == 2
        fsz = [fsz 1];
    end
    F = reshape(F,[fsz(1) prod(fsz(2:3))]);
    %F = bsxfun(@minus,F,mean(F,2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fprintf(['Starting fourier calculation.\n']);
    [M] = mfftm(N(2),NF);
    fF = mtimesx(M,F);
    fF = reshape(fF,[size(M,1) fsz(2:3)]);
    fF = permute(fF,[2 1 3]);
    flush = floor((0:42)*2*pi) + 2;
    flush(flush  > size(fF,2)) = [];
    %{
    for r = 1:numel(flush)
        fF(r,flush(r):end) = 0;
    end
    %}
    
    AMP = abs(fF);
    ANG = angle(fF);
    
    for k = 1:size(ANG,2)
        for r = 1:size(ANG,1)
            f = exp(1i*ANG(r,k));
        end
    end
    
    F = exp(1i*ANG);
    DIS = conj(F(1:end-1,:)).*F(2:end,:);
    
    %{
    ANG = unwrap(ANG,[],1);
    dA = diff(ANG,1,1);
    odA = -sign(dA).*(2*pi - abs(dA));
    dA = cat(3,dA,odA);
    dA = sort(dA,3,'descend');
    dA = dA .^ .5;
    %}
    
    
    AMP = mean(AMP,3);
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove this attempt
    aT = -pi/4;
    dA = diff(ANG,1,1);
    deltaTheta = dA;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % at choice theta
    choiceTheta = -sign(deltaTheta).*2*pi + deltaTheta;
    d = cat(3,deltaTheta,choiceTheta);
    d = sort(d,3);
    cur = ANG(1,:,1);
    whole = cur;
    
    target = 0:(size(d,1)-1);
    kvec = 0:(size(d,2)-1);
    kvec = ones(1,size(d,2));
    %target = target.^2;
    target = target'*kvec;
    
    for e = 1:(size(d,1))
        pot = bsxfun(@plus,cur,d(e,:,:));
        for c = 1:size(pot,2)
            [~,nidx(c)] = min(abs(pot(1,c,:)-target(e,c)),[],3);
            cur(1,c) = pot(1,c,nidx(c));
        end
        whole = [whole;cur];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    
    
    %{
    odA = 2*pi + dA;
    dA(dA < aT) = odA(dA < aT);
    dA = [ANG(1,:);dA];
    %}
    
    %dA = [ANG(1,:);whole];
    
    %ANG = cumsum(dA,1);
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nANG = ANG.*ikm;
    DIS = zeros(size(ANG,2),size(ANG,2),size(ANG,1),size(ANG,3));
    for m = 1:size(ANG,3)
        for rho = 1:size(ANG,1)
            tmpM = squeeze(nANG(rho,:,m))';
            tmpM = squareform(pdist(tmpM));
            DIS(:,:,rho,m) = 1*dK + 2*pi*abs(tmpM);
        end
    end
    DIS = permute(DIS,[3 1 2 4]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}

    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DIS = zeros(size(ANG,2),size(ANG,2),size(ANG,1),size(ANG,3));
    for m = 1:size(ANG,3)
        for rho = 1:size(ANG,1)
            tmpM = squeeze(ANG(rho,:,m))';
            tmpM = pdist(tmpM);
            %tmpM = squareform(pdist(tmpM));
            %{
            choiceTheta = -sign(deltaTheta).*2*pi + sign(deltaTheta).*deltaTheta;
            choiceVector = [deltaTheta;choiceTheta];
            [~,idx] = min(abs(choiceVector),[],1);
            DIS = choiceVector(idx);
            %}
            DIS(:,rho,m) = tmpM;
        end
    end
    %DIS = permute(DIS,[3 1 2 4]);
    DIS = permute(DIS,[2 1 3]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nAMP = AMP;
    nAMP(:,1) = mean(nAMP(:));
    nAMP = bindVec(nAMP);
    NOR = zeros(size(nAMP,2),size(nAMP,2),size(nAMP,1));
    df = diag(ones(size(nAMP,2),1));
    df = -(df - 1);
    for e = 1:size(nAMP,1)
        tmp = nAMP(e,:)'*nAMP(e,:);
        tmp(:,1) = 0;
        tmp(1,:) = 0;
        tmp = tmp.*df;
        NOR(:,:,e) = tmp;
    end
    NOR = permute(NOR,[3 1 2]);
    NOR = NOR.^.5;
    DIS = DIS.*NOR;
    %}
    %{
    iANG = ANG(:,:,1);
    for e = 1:size(ANG,3)
        ANG(:,:,e) = ANG(:,:,e) - iANG;
    end
    %}
    %ANG = unwrap(ANG,[],3);
    ANG = ANG(:,:,1);
    %ANG = unwrap(ANG,[],1);
    %[DIS] = calcDM(ANG);
    DIS = mean(DIS,4);
    Z = AMP.*exp(1i*ANG);
   

    
end