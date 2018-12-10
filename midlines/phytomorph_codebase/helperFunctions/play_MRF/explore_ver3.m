T = permn([0 1],9);
P = zeros(size(T,1),1);
for e = 1:size(T,1)
    s = reshape(T(e,:),[3 3]);
    m = imfilter(s,ones(2));
    m(:,end) = [];
    m(end,:) = [];
    
    if any(m(:) >= 3)
        P(e) = 0;
    else
        if T(e,5) == 0
            P(e) = 0;
        else
            if sum(T(e,:)) ~= 3
                P(e) = 0;
            else
                P(e) = 1;    
            end
            
        end
    end
    %{
    if (sum(s(:)) == 1) & (s(5) == 0)
        P(e) = 1;
    end
    %}
end

gb = [[1 1 0];[0 1 0];[1 1 0]];
for e = 1:4
    gb = imrotate(gb,90);
    gbIDX = (2*ones(1,9)).^fliplr((0:8))*gb(:)+1;
    P(gbIDX) = 1;
end


gb = reshape(T(27,:),[3 3]);
for e = 1:4
    gb = imrotate(gb,90);
    gbIDX = (2*ones(1,9)).^fliplr((0:8))*gb(:)+1;
    P(gbIDX) = 1;
end

%%
%b = vTest(ones(25),P);
%b = vTest(eye(3),P);
%I = mag(eye(3),3);
%I
[baseFilter f] = generateFilterSet_v2(T,P);
%% simple shape 1
tI = zeros(5);
tI(8) = 1;
tI(12) = 1;
tI(14) = 1;
tI(18) = 1;

%% circle
for r = 1%:1000
    %{
    tI = strel('disk',11,0);
    tI = padarray(getnhood(tI),[2 2],0,'both');
    tI = edge(tI);
    tI = bwmorph(tI,'skeleton',inf);
    
    vTest(double(tI),P);
    tI = double(tI);
    %}

    
    tI = zeros(5);
    tI(8) = 1;
    tI(12) = 1;
    tI(14) = 1;
    tI(18) = 1;
    
    tI(7) = 1;
    tI(9) = 1;
    tI(17) = 1;
    tI(19) = 1;
    
    
    %imshow(tI,[]);
    %waitforbuttonpress
    
    
    
    tar = zeros(5);
    for mg = 1:5
        tar = mag(tar,3);        
    end
    tar((end-1)/2,(end-1)/2) = 1;
    dtar = bwdist(tar);
    tar = dtar < 5*(3^4)/2;
    [g1 g2] = gradient(double(tar));
    tar = (g1.^2 + g2.^2).^.5;
    imshow(tar)
    for mg = 1:5
        tl{mg} = abs(imresize(tar,5*[3 3].^mg,'bilinear'));
        %tl{mg} = abs(imresize(tar,5*[3 3].^mg));
    end
    
    
    
    
    
    close all
    for l = 1:4
        kidx = find(tI);
        %kidx = kidx(randperm(numel(kidx)));
        [i1 i2] = find(tI);
        FM = ones(2^9,numel(kidx));
        IDX = zeros(size(tI));
        IDX = mag(IDX,3);
        IDX = reshape(1:numel(IDX),size(IDX));

        nI = zeros(size(IDX));

        for e = 1:numel(kidx)
            %e = 1;
            tmp = zeros(size(tI));
            tmp(kidx(e)) = 1;
            tmp = mag(tmp,3);
            tidx = find(tmp);
            sq = tI((i1(e)-1):(i1(e)+1),(i2(e)-1):(i2(e)+1));
            
            selectPoint = zeros(size(tI));
            selectPoint(kidx(e)) = 1;
            selectPoint = mag(selectPoint,3);
            selidx = find(selectPoint);
            ttTar = tl{l};
            tTar = [];
            for r = 1:4
                tTar(:,:,r) = imrotate(ttTar,(r-1)*90);
            end
            tTar = mean(tTar,3);
            reshape(tTar(selidx),[3 3])
            [of d new_image] = replace_v2(sq,FM(:,e),baseFilter,f,T,P,tTar,selidx);
            reshape(new_image,[3 3])
            nI(tidx) = new_image;

            
            for pl = 1:numel(of)
                nidx = sub2ind(size(tI),i1(e)+d(pl,1),i2(e)+d(pl,2));
                nidx = find(kidx==nidx);
                nidx
                e
                FM(:,nidx) = FM(:,nidx) & of{pl};
            end
            %{
            nidx = sub2ind(size(tI),i1(e)+d1(1),i2(e)+d1(2));
            nidx = find(kidx==nidx);
            FM(:,nidx) = FM(:,nidx) & of1;


            nidx = sub2ind(size(tI),i1(e)+d2(1),i2(e)+d2(2));
            nidx = find(kidx==nidx);
            FM(:,nidx) = FM(:,nidx) & of2;
            %}
            e;
            %imshow(nI,[])
            %if ~vTest(nI,P)
            %    break
            %end
            %waitforbuttonpress
            
            
            imshow(nI,[]);
            drawnow
            waitforbuttonpress
            
            
        end
        %nI = double(bwmorph(nI,'skel',inf));
        tI = nI;
        imshow(nI,[]);
        drawnow
        waitforbuttonpress
    end
    imshow(tI,[]);
    drawnow
    G{r} = tI;
    title(num2str(vTest(G{r},P,0)))
    waitforbuttonpress
    drawnow
end


%%
MTI = cat(2,logical(T),~logical(T));
fidx = find(P);
%%
for e = 1:numel(fidx)
    t = T(fidx(e),:);
    t = reshape(t,[3 3]);
    imshow(t,[]);
    t
    waitforbuttonpress
end

%%
pidx1 = find(sE(:)==1);
pidx0 = find(sE(:)==0);
v0 = E(pidx0);
v1 = E(pidx1);
[f0,xi0] = ksdensity(v0(1:10000),linspace(0,50,500));
[f1,xi1] = ksdensity(v1(1:10000),linspace(0,50,500));
%%
midx0 = T(:,5) == 0;
midx1 = T(:,5) == 1;
%%
tmp = E(:,:,1);
tmp = imresize(tmp,.15);
tmp = imcrop(tmp,[]);

%{
tmp = [[0 1 0];[1 0 1];[0 1 0]];
tmp = padarray(tmp, [30 20], 0, 'both');
%}

SZ = size(tmp);
tmp = padarray(tmp, [1 1], 'replicate', 'both');


p0 = interp1(xi0,f0,tmp,'linear',eps);
p1 = interp1(xi1,f1,tmp,'linear',eps);



%{
p1 = .8*tmp;
%p0 = double(~logical(p1));
p0 = imcomplement(p1);
%}

T = p0 + p1;
p1 = p1 .* T.^-1;
p0 = p0 .* T.^-1;


p_O_0 = p0;
p_O_1 = p1;


%p_O_1 = imcomplement(p_O_0);


midx0 = P == 0;
midx1 = P == 1;
close all
PS = [];



recent = 0;


for e = 1:2000000
    
    
    [snp0 snp1] = update(p0,p1,MTI,midx0,midx1,SZ,p_O_0,p_O_1);
    
    
    
    imshow(p1,[]);
    drawnow
   
    
    e
    
    
    
    deltaP = abs(p1 - snp1);
    
    pidx(e) = recent(1);
    while sum(recent == pidx(e)) ~= 0
        [V(e) pidx(e)] = max(deltaP(:));
        deltaP(pidx(e)) = 0;
    end
    
    
    p1(pidx(e)) = snp1(pidx(e));
    p0(pidx(e)) = snp0(pidx(e));
    %{
    MD  = (snp1 - p1);
    
    
    
    c1 = p1;
    c0 = p0;
    Z = zeros(size(p1));
    Z(pidx(e)) = 1;
    
    fidx = 1;
    while ~isempty(fidx)
        %z = 1:100
        tZ = imdilate(Z,strel('disk',z,0));
        
        [tsnp0 tsnp1] = update(c0,c1,MTI,midx0,midx1,SZ,p_O_0,p_O_1);
        
        mD = (tsnp1 - c1);
        
        fD = (mD - MD);
        
        
        fidx = find(tZ.*fD ~= 0);
        
        
        %c0(fidx) = tsnp0(fidx);
        %c1(fidx) = tsnp1(fidx);
        
        c0(fidx) = c0(fidx) - fD(fidx);
        c1(fidx) = c1(fidx) + fD(fidx);
        if ~isempty(fidx)
            imshow(c1,[]);
            %waitforbuttonpress
            drawnow;
        end
    end
    
    %}
    
    Z = zeros(size(p1));
    Z(pidx(e)) = 1;
    Z = imdilate(Z,ones(11));
    out = flattenMaskOverlay(p1,logical(Z),1);
    imshow(out,[]);
    drawnow
    %waitforbuttonpress
    
    recent = [recent pidx(e)];
    if numel(recent) > 100
        recent(1) = [];
    end
    pidx(e)
end
%%
close all
for e = 1:size(PS,3)
    imshow(PS(:,:,e),[]);
    drawnow
    waitforbuttonpress
end
%%
close all
for e = 1:size(PS,3)
    imshow(PS(:,:,e),[])
    drawnow
    pause(.5)
end




























