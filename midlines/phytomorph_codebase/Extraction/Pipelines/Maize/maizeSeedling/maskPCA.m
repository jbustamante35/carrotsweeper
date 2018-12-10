clear all
FilePath = '/home/nate/Downloads/imageMasks/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
%

parfor e = 1:numel(FileList)
    tic
    R{e} = getSingleMask(FileList{e});
    toc
end
%%
kidx = kmeans(stack(1:50,:)',2);



%% profile
stack = zeros([2000 Tcnt]);
cnt = 1;
close all
rIDX = [];
for r = 1:numel(R)
    for e = 1:numel(R{r})
        profile = flip(sum(R{r}(e).Image,2),1);
        tmp = zeros(2000,1);
        tmp(1:numel(profile)) = profile;
        stack(:,cnt) = tmp;
        cnt = cnt + 1;
        %stack = [stack tmp];
        rIDX = [rIDX;[r e]];
        %{
        stack(1:numel(profile)) = stack(1:numel(profile)) + profile;
        cnt = cnt + 1;
        plot(stack*cnt^-1);
        drawnow
        %}
        %errorbar(mean(stack,2),std(stack,1,2)*size(stack,2)^-.5);
        %drawnow
    end
    r
end
%% count number of images
Tcnt = 1;
mSZ = [0 0];
for r = 1:numel(R)
    for e = 1:numel(R{r})
        Tcnt = Tcnt + 1;
        SZ = size(R{r}(e).Image);
        mSZ(1) = max(SZ(1),mSZ(1));
        mSZ(2) = max(SZ(2),mSZ(2));
    end
end
%%
S = zeros(500,750,Tcnt,'single');
%
cnt = 1;
for r = 1:numel(R)
    for e = 1:numel(R{r})
        tic
        [S(:,:,cnt) w(cnt)] = centerPlant(R{r}(e).Image);
        cnt = cnt + 1;
        toc
        
        
    end
end
%{
%%
SBK = S;
%%
%}
S(:,:,w==0) = [];
%
rS = zeros(250,375,size(S,3),'single');
for e = 1:size(S,3)
    tic
    rS(:,:,e) = imresize(S(:,:,e),[250 375]);
    toc
end
%%
sz = size(rS);
%tS = reshape(cat(3,rS,flipdim(rS,2)),[prod(sz(1:2)) 2*sz(3)]);
tS = reshape(rS,[prod(sz(1:2)) sz(3)]);
%
close all
[fS fC fU fE fL fERR fLAM] = PCA_FIT_FULL_T(double(tS),10);
U = reshape(fU,sz(1:2));
imshow(U,[]);
%% sweep
close all
for d = 1:3
    uu = mean(fC,2);
    R = std(fC(d,:),1,2);
    l = linspace(-R,R,20);
    for sw = 1:numel(l)
        tmp = uu;
        tmp(d) = l(sw);
        M = PCA_BKPROJ_T(tmp,fE,fU);
        M = reshape(M,sz(1:2));
        imshow(M,[]);
        title(num2str(d));
        drawnow
        pause(.1);
    end
end
%%
close all
L1 = linspace(min(fC(1,:)),max(fC(1,:)),20);
clear uC
grandU = [];
for e1 = 1:(numel(L1)-1)
    %plot3(fC(1,:),fC(2,:),fC(3,:),'.')
    hold on
    idx = find(fC(1,:) > L1(e1) & fC(1,:) < L1(e1+1));
    if e1 == 5
        plot3(fC(1,idx),fC(2,idx),fC(3,idx),'r.');
    end
    uC(:,e1) = mean(fC(:,idx),2);
    
    [vS vC{e1} vU{e1} vE{e1} vL vERR vLAM] = PCA_FIT_FULL_T(double(fC(1:3,idx)),3);
    
    
    grandU(:,e1) = mean(tS(:,idx),2);
    IDX{e1} = idx;
    
    
    plot3(uC(1,e1),uC(2,e1),uC(3,e1),'k.');
    %hold off
    drawnow
end



plot3(uC(1,:),uC(2,:),uC(3,:),'k')

dU = diff(uC,2,1);



[TAN] = gradient(uC(1:3,:));
[NOR] = gradient(TAN(1:3,:));

for e = 1:size(dU,2)
    TAN(:,e) = TAN(:,e) / norm(TAN(:,e));
    
    NOR(:,e) = NOR(:,e) - NOR(:,e)'*TAN(:,e)*TAN(:,e);
    NOR(:,e) = NOR(:,e) / norm(NOR(:,e));
    dU(:,e) = dU(:,e) / norm(dU(:,e));
    
    
    
    
    for k = 1:3
        if vE{end}(:,k)'*vE{e}(:,k) < 0
            vE{e}(:,k) = -vE{e}(:,k);
        end
    end
   
   
    hold on
    quiver3(uC(1,e),uC(2,e),uC(3,e),vE{e}(1,1),vE{e}(2,1),vE{e}(3,1),'r');
    quiver3(uC(1,e),uC(2,e),uC(3,e),vE{e}(1,2),vE{e}(2,2),vE{e}(3,2),'g');
    quiver3(uC(1,e),uC(2,e),uC(3,e),vE{e}(1,3),vE{e}(2,3),vE{e}(3,3),'b');
    
end

for e = 2:size(dU,2)
    if NOR(:,e-1)'*NOR(:,e) < 0
        NOR(:,e) = -NOR(:,e);
    end
    
     quiver3(uC(1,e),uC(2,e),uC(3,e),TAN(1,e),TAN(2,e),TAN(3,e),'m');
    quiver3(uC(1,e),uC(2,e),uC(3,e),NOR(1,e),NOR(2,e),NOR(3,e),'c');
end

%%
close all
figure
exD = 1;
mag = 1;
MS = [];
for e = 1:size(grandU,2)
    gM = [];
    tmpM = reshape(grandU(:,e),sz(1:2));
    gM = [gM tmpM];
    
    
    [tmpC] = PCA_REPROJ_T(double(fC(1:3,IDX{e})),vE{e},vU{e});
    J = mean(tmpC,2);
    sJ = std(tmpC(exD,:),1,2);
    LN = linspace(-10*sJ,10*sJ,3);
    
    [tmpC] = PCA_REPROJ_T(double(fC(1:3,IDX{e})),NOR(:,e),vU{e});
    J = mean(tmpC,2);
    sJ = std(tmpC(exD,:),1,2);
    LN = linspace(-mag*sJ,mag*sJ,3);
    
    for l = 1:numel(LN)
        tmpJ = J;
        tmpJ(exD) = LN(l);
        tmpJ = PCA_BKPROJ_T(tmpJ,NOR(:,e),vU{e});
        %tmpJ = PCA_BKPROJ_T(tmpJ,vE{e},vU{e});
        tmpy = mean(fC,2);
        tmpy(1:3)
        tmpy(1:3) = tmpJ;
        tmpy(1:3)
        tmpy = PCA_BKPROJ_T(tmpy,fE,fU);
        tmpy = reshape(tmpy,sz(1:2));
        gM= [gM tmpy];
       
    end
    
    imshow(gM,[]);
    drawnow
    waitforbuttonpress
    MS = [MS;gM];
    imshow(tmpM,[])
    drawnow
    pause(.1)
end
close all
imshow(MS,[]);
%% skeleton 
for e = 1:size(S,3)
    tic
    [LP{e}] = massiveP(S(:,:,e));
    toc
end
%%
for e = 1:numel(LP)
    LP{e}(1,1)
end
%%
for e = 1:numel(LP)
    len(e) = size(LP{e},1);
end
DIS = zeros(8,8,8);
LT = [[-1 0];[-1 1];[0 1];[1 1];[1 0];[1 -1];[0 -1];[-1 -1]];
for e = 1:numel(LP)
    try
        d = diff(LP{e});
        for p = 1:size(d,1)
            f = find(all(bsxfun(@eq,LT,d(p,:)),2));
            CC{e}(p) = f;
        end
    catch
        CC{e} = [];
    end
end
%% 
for e = 1:numel(CC)
    tic
    for l = 1:(numel(CC{e})-2)
        DIS(CC{e}(l),CC{e}(l+1),CC{e}(l+2)) = DIS(CC{e}(l),CC{e}(l+1),CC{e}(l+2)) + 1;
    end
    toc
end
%%
len = [];
for e = 1:numel(CC)
    len(e) = numel(CC{e});
end
[ss,ridx] = sort(len);
CC(len < 30) = [];
W = 10;
for l = min(len):max(len)
    fidx = find((len >= (l - W)) &  (len <= (l + W)));
    for f = 1:numel(fidx)
    
    end
end
%%

Z = zeros(size(S,1),size(S,2));
chain = [3 3];
h = 1;
for l = 1:300
    slab = squeeze(DIS(chain(h),chain(h+1),:));
    slab = slab / sum(slab);
    cs = cumsum(slab);
    r =  rand(1);
    fidx = find(cs < r);
    r = fidx(end) + 1;
    chain = [chain r];
end
dP = [375,500];
for c = 1:numel(chain)
    dP = [dP ; LT(chain(c),:).*[1 -1]];
end
P = cumsum(dP,1);
for e = 1:size(P,1)
    Z(P(e,2),P(e,1)) = 1;
end
close all
plot(P(:,2),P(:,1))
imshow(Z,[])





