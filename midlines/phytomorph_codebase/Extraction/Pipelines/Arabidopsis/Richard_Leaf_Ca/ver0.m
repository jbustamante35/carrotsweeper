%% read data
javaaddpath('/mnt/scratch1/phytomorph_dev/helperFunctions/JARs/loci_tools.jar')
fileList = '/home/nate/Downloads/r.nd2';
imreadBFmeta(fileList)
D = imreadBF(fileList,1,[1:788],1);
E = imreadBF(fileList,1,[1:788],2);
%% fix zero
for e = 1:size(D,3)
    fidx = find(D(:,:,e)==0);
    tmp = D(:,:,e);
    tmp(fidx) = mean(tmp(:));
    D(:,:,e) = tmp;
    e
end
%% make ratio
close all 
R = zeros(size(E));
for e = 1:size(D,3)
    R(:,:,e) = E(:,:,e).*D(:,:,e).^-1;
    e
end
%% perform crop
[J,BOX] = imcrop(R(:,:,1),[]);
%% filter stack
close all
J = zeros([size(J) size(R,3)]);
for e = 1:size(R,3)
    J(:,:,e) = imcrop(R(:,:,e),BOX);
    J(:,:,e) = imfilter(J(:,:,e),fspecial('disk',5),'replicate');
    e
end
%% dont need
close all
for e = 1:size(J,3)
    imshow(J(:,:,e),[1 3]);
    drawnow
end
%% dont need
close all
S = std(J,1,3);
U = mean(J,3);
imshow(U,[])
%% dont need
[x y v] = impixel(U,[]);
d = [];
for e = 1:numel(x)
    d(:,e) = squeeze(J(y(e),x(e),:));
end
close all
plot(d)
%% dont need
fr = S;
[P sidx] = max(fr(:));
Z = zeros(size(fr));
Z(sidx) = 1;
Z = imdilate(Z,strel('disk',21));
imshow(Z,[])
%% dont need
close all
imshow(max(J,[],3),[])
%%  dont need
sig = reshape(J,[size(J,1)*size(J,2) size(J,3)]);
%% dont need
plot(sig(1:1000:end,:)');
%% smooth in time
sig = reshape(J,[size(J,1)*size(J,2) size(J,3)]);
sig = imfilter(sig,fspecial('average',[1 11]),'replicate');
sJ = reshape(sig,size(J));
%% derivative NEED
dJ = diff(sJ,1,3);
%[d1 d2 d3] = gradient(J);
%dJ = (d1.^2 + d2.^2 + d3.^2).^.5;
%% try smoothing dJ
sdJ = zeros(size(sdJ));
for e = 1:size(dJ,3)
    sdJ(:,:,e) = imfilter(dJ(:,:,e),fspecial('disk',21),'replicate');
    e
end
%% for showing derivative
close all
tmpS = [];
for e = 1:size(dJ,3)
    tmpS(:,:,e) = imresize(sdJ(:,:,e),.5);
    e
end
%% show derivative
close all
for e = 1:size(tmpS,3)
    imshow(tmpS(:,:,e),[])
    drawnow
end
%%
[d1 d2 d3] = gradient(sdJ);
%tough = (d1.^2 + d2.^2 + d3.^2).^.5;
tough = (d1.^2 + d2.^2);
%% second der for wave
[~,~,DD3] = gradient
%% 
L = fspecial('laplacian');
for t = 1:size(J,3)
    LS(:,:,t) = imfilter(sJ(:,:,t),L,'replicate');
end
%%
vel = DD3.*LS(:,:,1:end-1).^-1;
vel = vel.^.5;

%% show local max
close all
for e = 1:size(G,3)
    imshow(real(LS(:,:,e)),[])
    drawnow
end
%% show tough
close all
for e = 1:size(tough,3)
    imshow(tough(:,:,e),[])
    drawnow
end
%% original 
close all
for e = 1:size(J,3)
    imshow(sdJ(:,:,e),[])
    drawnow
end
%% show K
close all
for e = 1:size(K,4)
    tmp = prod(K(:,:,:,e),3);
    imshow(tmp,[])
    drawnow
end
%% show der
close all
for e = 1:size(sdJ,3)
    MSK = sdJ(:,:,e) > .015;
    imshow(sdJ(:,:,e).*MSK,[])
    %imshow(sdJ(:,:,e).*MSK,[])
    drawnow
end
%%
close all
for l = 1:10
    for e = 1:size(J,3)
        MSK = sdJ(:,:,e) > .013;
        dB = bwboundaries(MSK);
        %out = flattenMaskOverlay((J(:,:,e)-M3)/M4,MSK);

        imshow((J(:,:,e)-M3)/M4,[])
        hold on
        for i = 1:numel(dB)
            plot(dB{i}(:,2),dB{i}(:,1),'r')
        end
        hold off
        drawnow
    end
end
%%
BLOB = sdJ > .013;
CC = bwconncomp(BLOB);
%% explore single space-time event
close all
figure
hold on
for k = 1:CC.NumObjects
    [w1 w2 w3] = ind2sub(size(BLOB),CC.PixelIdxList{k});
    plot3(w1,w2,w3,'.')
    eSZ(k) = numel(w1);
end
%% look at largest event
close all
[~,midx] = max(eSZ);
[w1 w2 w3] = ind2sub(size(BLOB),CC.PixelIdxList{midx});
plot3(w1,w2,w3,'.')
%%
singleBLOB = zeros(size(BLOB));
singleBLOB(CC.PixelIdxList{midx}) = 1;
F = bwdist(singleBLOB);
G = bwdist(~singleBLOB);
H = F - G;
%%
Hs = smooth3(H,'gaussian',[9 9 9]);
%% make surface ver0 of largest event
fv = isosurface(H,0);
[g1 g2 g3] = gradient(H);
GM = (g1.^2 + g2.^2 + g3.^2).^-.5;
gn1 = g1.*GM;
gn2 = g2.*GM;
gn3 = g3.*GM;
%% 
for e = 1:size(H,3)
    [n1(:,:,e) n2(:,:,e)] = gradient(H(:,:,e));
    e
end
N = (n1.^2 + n2.^2).^-.5;
nn1 = n1.*N;
nn2 = n2.*N;
t1 = nn2;
t2 = -nn1;

nn1 = double(nn1);
nn2 = double(nn2);
%%
p = patch(fv);
set(p,'facecolor','red','edgecolor','none');
%% movie with largest event
close all
h1 = figure;
h2 = figure;
SP = [];
for e= 1:size(J,3)
    figure(h1);
    imshow(J(:,:,e),[1 4]);
    Z = zeros(size(J,1),size(J,2));
    Zn = zeros(size(J,1),size(J,2));
    Zl = zeros(size(J,1),size(J,2));
    hold on
    fidx = find(w3==e);
    
    
    if ~isempty(fidx)
        nidx = find(w3==(e+1));
        lidx = find(w3==(e-1));
        ind = sub2ind(size(Z),w1(fidx),w2(fidx));
        nind = sub2ind(size(Z),w1(nidx),w2(nidx));
        lind = sub2ind(size(Z),w1(lidx),w2(lidx));
        Z(ind) = 1;
        Zn(nind) = 1;
        Zl(lind) = 1;
        dB = bwboundaries(Z);
        dBn = bwboundaries(Zn);
        dBl = bwboundaries(Zl);
        hold on
        
        for k = 1:numel(dBn)
            plot(dBn{k}(:,2),dBn{k}(:,1),'w');
        end
        
        
        for k = 1:numel(dBl)
            plot(dBl{k}(:,2),dBl{k}(:,1),'y');
        end
        
        
        
        
        
        for k = 1:numel(dB)
            plot(dB{k}(:,2),dB{k}(:,1),'r');
            
             
            v1 = ba_interp2(n1(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),'cubic');
            v2 = ba_interp2(n2(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),'cubic');
            %{
            V = (v1.^2 + v2.^2).^-.5;
            v1 = v1.*V;
            v2 = v2.*V;
            v3 = zeros(size(v1));
            %}
            W1 = ba_interp2(g1(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),'cubic');
            W2 = ba_interp2(g2(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),'cubic');
            W3 = ba_interp2(g3(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),'cubic');
            %{
            W = (W1.^2 + W2.^2 + W3.^2).^-.5;
            W1 = W1.*W;
            W2 = W2.*W;
            W3 = W3.*W;
            %}
            
            t1 = v2;
            t2 = -v1;
            t3 = v3;
            
            A = [t1 t2 t3];
            B = [W1 W2 W3];
            C = cross(A,B,2);
            speed = bsxfun(@times,C,C(:,3).^-1);
            
            quiver(dB{k}(:,2),dB{k}(:,1),speed(:,1),speed(:,2),'c')
            speed = sum(speed(:,1:2).*speed(:,1:2),2).^.5;
            %speed = dot(speed(:,1:2),speed(:,1:2),2);
            
            cp = mean(dB{k},1);
           
            
            
            SP = [SP;speed];
            
            
            
            quiver(dB{k}(:,2),dB{k}(:,1),v1,v2,'g');
            text(cp(2),cp(1),num2str(mean(speed),'% 10.2f'),'backgroundcolor','w');
            
            %quiver(dB{k}(:,2),dB{k}(:,1),t1,t2,'b');
        end
        
    end
    %figure(h2)
    %ksdensity(speed(:));
    %plot(w2(fidx),w1(fidx),'.')
    hold off
    drawnow
end
%%
para.scales.value = [5];
para.resize.value = 1;
U = mean(J,3);
[u1 u2] = gradient(U);
dU = (u1.^2 + u2.^2).^.5;
dU = surKur(U,para);
v = imcomplement(dU(:,:,2));
imshow(v,[]);
%% show composite

M1 = min(sdJ(:));
M2 = max(sdJ(:)-M1);
M3 = min(J(:));
M4 = max(J(:)-M3);
M5 = min(LS(:));
M6 = max(LS(:) - M5);
%%
vn = bindVec(v);
close all
for e = 1:size(sdJ,3)
    tmp1 = (sdJ(:,:,e) - M1)/M2;
    tmp2 = (J(:,:,e) - M3)/(M4*.9);
    %vn = (LS(:,:,e) - M5)/M6;
    tmp = cat(3,tmp1,vn,tmp2);
    COMP(:,:,:,e) = tmp;
    imshow(tmp,[])
    drawnow
end
%% for 
for r = 1:size(COMP,4)
    RE(:,:,:,r) = imresize(COMP(:,:,:,r),.6);
    e
end
%% f
close all
for loop = 1:10
    for t = 1:2:size(COMP,4)
        
        imshow(RE(:,:,:,t),[]);
        title(num2str(t));
        drawnow
    end
end
%% K
para.scales.value = [5];
para.resize.value = .15;
K = zeros([size(sdJ,1) size(sdJ,2) 2 size(sdJ,3)]);
for e = 1:size(sdJ)
    K(:,:,:,e) = surKur(sdJ(:,:,e),para);
    e
end

%% find first puff
close all
DELTA = 50;
[K,idx] = max(sdJ(:,:,DELTA:end),[],3);
%% find max for all puffs
close all
SZ = 71;
TMS = 5;
[mG] = minmaxfilt(sdJ,[SZ SZ TMS], 'max', 'same');
G = (mG == sdJ);
gidxM = find(G);
[yp xp TM] = ind2sub(size(G),gidxM);
%% for finding first
close all
[yp xp] = find(imerode(idx,strel('disk',21)) == idx);
lidx = find(imerode(idx,strel('disk',21)) == idx);
TM = idx(lidx);
%[yp xp] = find(imdilate(K,strel('disk',101)) == K);
Z = zeros(size(idx));
for e = 1:numel(xp)
    Z(yp(e),xp(e)) = 1;
end
imshow(imdilate(Z,strel('disk',10)),[]);
%% filter on border
Z = zeros(size(dJ,1),size(dJ,2));
for r = 1:4
    Z(1:SZ,:) = 1;
    Z = imrotate(Z,90);
end
imshow(Z,[])
sam = [];
for e = 1:numel(xp)
    sam(e) = Z(yp(e),xp(e));
end
rm = find(sam==1);
xp(rm) = [];
yp(rm) =[];
TM(rm) = [];
%%
BKxp = xp;
BKyp = yp;
BKTM = TM;
%%
xp = BKxp;
yp = BKyp;
TM = BKTM;
%% sample events
samp = zeros(size(TM));
for e = 1:numel(xp)
    gidx = find((TM==e));
    if ~isempty(gidx)
        samp(gidx) = interp2(mG(:,:,e),xp(gidx),yp(gidx));
    end
end
T = graythresh(samp);
rm = samp < T;
xp(rm) = [];
yp(rm) = [];
TM(rm) = [];
%%
close all
for e = 1:size(tmpS,3)
    %imshow(sdJ(:,:,e),[]);
    imshow(J(:,:,e),[1 3]);
    gidx = find((TM+TMS==e));
    hold on
    %plot(xp,yp,'b*')
    plot(xp(gidx),yp(gidx),'r*')
    hold off
    drawnow
    %waitforbuttonpress
end
%%
close all
for e = 1:size(mG,3)
    tmp1 = bindVec(mG(:,:,e));
    tmp2 = bindVec(J(:,:,e));
    tmp3 = cat(2,tmp1,tmp2);
    imshow(tmp3,[])
    drawnow
end
