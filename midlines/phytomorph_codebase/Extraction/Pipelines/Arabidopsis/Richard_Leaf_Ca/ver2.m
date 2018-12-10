FilePath = '/mnt/snapper/nate/forRichard/forACA/';
FileList = {};
FileExt = {'nd2'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% add jar
javaaddpath('/mnt/scratch1/phytomorph_dev/helperFunctions/JARs/loci_tools.jar')
%% read e

for e = 1:numel(FileList)
    meta = imreadBFmeta(FileList{e});
    MF = meta.nframes;
    qM(e).TOT = MF;
    for c = 1:2
        
        qM(e).D(:,:,1,c) = imreadBF(FileList{e},1,[1],c);
        qM(e).D(:,:,2,c) = imreadBF(FileList{e},1,[round(MF/2)],c);
        qM(e).D(:,:,3,c) = imreadBF(FileList{e},1,[MF],c);
    end
    e
end
%% view data
close all
for e = 1:numel(qM)
    imshow(qM(e).D(:,:,:,1)/255,[])
    drawnow
    pause(.3)
end
%%
%% view data
close all
for e = 1:numel(qM)
    d1 = qM(e).D(:,:,3,1) - qM(e).D(:,:,1,1);
    d2 = qM(e).D(:,:,3,2) - qM(e).D(:,:,1,2);
    dU(e,:) = [mean(d1(:)) mean(d2(:))];
end
%%
close all
kidx = kmeans(abs(dU(:,1)),3);

options = statset('Display','iter');
GMModel = fitgmdist(dU(:,1),3,'Options',options,'RegularizationValue',0.01);
kidx = GMModel.cluster(dU(:,1));
for k = 3
    fidx= find(kidx==k);
    for f = 1:numel(fidx)
        imshow(qM(fidx(f)).D(:,:,:,1)/255,[]);
        title(num2str(k))
        drawnow
    end
end
%%
imshow(D(:,:,[1 50 100])/255,[])

%% read data
javaaddpath('/mnt/scratch1/phytomorph_dev/helperFunctions/JARs/loci_tools.jar')
fileList = '/home/nate/Downloads/r.nd2';
%fileList = '/mnt/spaldingdata/nate/160529_aca4.1-11.1 T1 plantE.nd2';
MF = 788;
meta = imreadBFmeta(fileList);
D = imreadBF(fileList,1,[1:MF],1);
E = imreadBF(fileList,1,[1:MF],2);
%%
fileList = '/home/nate/Downloads/r.nd2';
fileList = '/mnt/spaldingdata/nate/160529_aca4.1-11.1 T1 plantE.nd2';
MF = 100;
imreadBFmeta(fileList)
DL = imreadBF(fileList,1,[MF],1);
EL = imreadBF(fileList,1,[MF],2);
%% read leasion
fileList2 = '/mnt/spaldingdata/nate/160601_aca4.1-11.1 +65 T1 plantC.nd2';
imreadBFmeta(fileList)
D1 = imreadBF(fileList2,1,[100],1);
E1 = imreadBF(fileList2,1,[100],2);
R1(:,:,e) = E1(:,:,e).*D1(:,:,e).^-1;
%% register
[aerial_points,ortho_points] = cpselect(bindVec(R1(:,:,1)),bindVec(sJ(:,:,1)),'Wait',true);
t_concord = cp2tform(aerial_points,ortho_points,'projective');
RR = imtransform(bindVec(R1(:,:,1)), t_concord,'XData',[1 size(sJ,2)],'YData',[1 size(sJ,1)]);
%%
[cR BOX] = imcrop(RR,[]);
cJ = imcrop(sJ(:,:,1),BOX);
imshow(cat(3,bindVec(cR),zeros(size(cJ)),bindVec(cJ)));
%%
pos = {};
for e = 1:11
     imshow(cR,[]);
    for k = 1:numel(pos)
        plot(pos{k}(:,1),pos{k}(:,2),'r')
    end
    hold on
   
    h = imfreehand;
    pos{e} = wait(h);
end
%% crop out leaf
imshow(cJ,[]);
h = imfreehand;
leafR = wait(h);

leafM = poly2mask(leafR(:,1),leafR(:,2),size(cJ,1),size(cJ,2));
leafM = repmat(leafM,[1 1 3]);
%% overlay trace
figure
%cJ = imcrop(sJ(:,:,end),BOX);
%cJ = imcrop(Ca_map,BOX);
%cJ = imcrop(max(ssdJ,[],3),BOX);
%cJ = imcrop(sP,BOX);
cJ = imcrop(H,BOX);
%cJ = imcrop(P,BOX);
imshow(cJ.*leafM,[]);

%imshow(max(ssdJ,[],3),[]);
hold on
for k = 1:numel(pos)
    plot(pos{k}(:,1),pos{k}(:,2),'c')
end
hold on
%% 
close all
G = bindVec(max(ssdJ,[],3));
G = imfilter(G,fspecial('disk',51));
%H = cat(3,bindVec(sP),zeros(size(sP)),bindVec(max(ssdJ,[],3)));
H = cat(3,bindVec(sP),bindVec(v),G);
imshow(H,[]);
%%
rT = [];
for t = 1:size(sJ,3)
    rT(:,:,t) = imresize(sJ(:,:,t),.25);
    t
end
sig = reshape(rT,[size(rT,1)*size(rT,2) size(rT,3)]);
p = [];
sig = sig(:,1:4:end);
for e = 1:size(sig,1)
    p(e,:) = polyfit(1:size(sig,2),sig(e,:),1);
    e
    size(sig,1)
end
P = reshape(p(:,1),[size(rT,1) size(rT,2)]);
P = imresize(P,[size(sJ,1) size(sJ,2)]);
sP = imfilter(P,fspecial('disk',51),'replicate');
%% next pos
close all
for e = 1:numel(pos)
    MSK = poly2mask(pos{e}(:,1),pos{e}(:,2),size(cR,1),size(cR,2));
    midx = find(MSK);
    VEC(e,:) = [mean(cR(midx)) mean(cJ(midx))];
end
plot(VEC(:,1),VEC(:,2),'.')
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
%% view Ratio movie
for e = 1:size(R,3)
    imshow(R(:,:,e),[])
    drawnow
end
%% perform crop
close all
[J,BOX] = imcrop(R(:,:,1),[]);
%% crop and filter the images
close all
J = zeros([size(J) size(R,3)]);
for e = 1:size(R,3)
    J(:,:,e) = imcrop(R(:,:,e),BOX);
    J(:,:,e) = imfilter(J(:,:,e),fspecial('disk',5),'replicate');
    e
end
%% view cropped smooth Ratio movie
for e = 1:size(J,3)
    imshow(J(:,:,e),[])
    drawnow
end
%% smooth in time
sJ = timeSmooth(J,11);
dJ = DGradient(sJ,1,3);
sdJ = timeSmooth(dJ,11);
d2J = DGradient(sdJ,1,3);
sd2J = timeSmooth(d2J,11);
%% smoothing dJ
ssdJ = zeros(size(sdJ));
for e = 1:size(dJ,3)
    ssdJ(:,:,e) = imfilter(sdJ(:,:,e),fspecial('disk',21),'replicate');
    e
end
%% smoothing d2J
ssd2J = zeros(size(sd2J));
for e = 1:size(dJ,3)
    ssd2J(:,:,e) = imfilter(sd2J(:,:,e),fspecial('disk',21),'replicate');
    e
end
%% derivative NEED
%dJ = diff(sJ,1,3);


%[~,~,dJ] = gradient(sJ);
%[d1 d2 d3] = gradient(J);
%dJ = (d1.^2 + d2.^2 + d3.^2).^.5;
%% explore image
close all
[c r v] = impixel(sJ(:,:,1),[]);
eidx = sub2ind([size(sJ,1) size(sJ,2)],r,c);
close all
xscale = 1;
%plot((1:788)*xscale,bindVec(squeeze(sJ(r,c,:))),'k');hold on
plot((1:788)*xscale,bindVec(squeeze(ssd2J(r,c,:))),'k');hold on
plot((1:788)*xscale,bindVec(squeeze(ssdJ(r,c,:))),'r');hold on
plot((1:788)*xscale,(squeeze(sJ(r,c,:))),'g');hold on
%plot((1:788)*xscale,.5*bindVec(squeeze(BLOB(r,c,:))),'r');hold on
%plot((1:788)*xscale,.5*bindVec(squeeze(front(r,c,:))),'b');hold on
%axis([0 788*15/60/60 0 1]);
xlabel('Time (hr)')
ylabel('Image Intensity')
figure;
plot((1:788)*xscale,(squeeze(ssd2J(r,c,:))),'k');hold on
%plot(bindVec(squeeze(sdJ(r,c,:))),'b');
%plot(bindVec(squeeze(sd2J(r,c,:))),'g');
%% find peaks for first der
sig = reshape(sJ,[size(sJ,1)*size(sJ,2) size(sJ,3)]);
dsig = imdilate(sig,ones(1,61));
front0 = dsig == sig;
front0 = reshape(front0,size(sdJ));
%% find peaks for first der
sig = reshape(ssdJ,[size(sdJ,1)*size(sdJ,2) size(sdJ,3)]);
dsig = imdilate(sig,ones(1,61));
front = dsig == sig;
front = reshape(front,size(sdJ));
%% find peaks for second der
sig = reshape(ssd2J,[size(sdJ,1)*size(sdJ,2) size(sdJ,3)]);
dsig = imdilate(sig,ones(1,61));
front2 = dsig == sig;
front2 = reshape(front2,size(sdJ));
%% show front1 and front2
close all
MOV = zeros(size(sJ,1),size(sJ,2),3,size(sJ,4));

for e = 1:size(front2,3)
    tmp = sJ(:,:,e);
    out = flattenMaskOverlay(bindVec(tmp),logical(wfront(:,:,e)),.5,'b');
    out = flattenMaskOverlay(out,logical(wfront2(:,:,e)),.5,'r');
    out = flattenMaskOverlay(out,logical(wfront0(:,:,e)),.5,'g');
    imshow(out,[]);
    title(num2str(e))
    drawnow
    MOV(:,:,:,e) = out;
end
%% look at average time signal
sig = reshape(wfront2,[size(sJ,1)*size(sJ,2) size(sJ,3)]);
sig2 = reshape(sJ,[size(sJ,1)*size(sJ,2) size(sJ,3)]);
fidx = find(sig);
[f1 f2] = ind2sub(size(sig),fidx);
rm = f2 < 150;
f2(rm) = [];
f1(rm) = [];
fidx(rm) = [];
delta = f2 + 150;
rm = delta > size(sig,2);
delta(rm) = [];
f1(rm) = [];
f2(rm) = [];
fidx(rm) = [];
M = zeros(numel(f1),150);
for e = 1:numel(f1)
    M(e,:) = sig2(f1(e),(f2(e):f2(e)+149));
    e
end
nM = bsxfun(@minus,M,M(:,1));
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(nM,3);
wv1 = zeros(size(sig));
wv2 = zeros(size(sig));
wv3 = zeros(size(sig));
for e = 1:numel(f1)
    wv1(fidx(e)) = pC(e,1);
    wv2(fidx(e)) = pC(e,2);
    wv3(fidx(e)) = pC(e,3);
end
wv1 = reshape(wv1,size(sJ));
wv2 = reshape(wv2,size(sJ));
wv3 = reshape(wv3,size(sJ));

wv1 = wv1 - min(wv1(:));
wv1 = wv1 / max(wv1(:));
wv2 = wv2 - min(wv2(:));
wv2 = wv2 / max(wv2(:));
wv3 = wv3 - min(wv3(:));
wv3 = wv3 / max(wv3(:));
%% normalize to first event
nM = bsxfun(@minus,M,M(:,1));
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(nM,3);
[icasig, A, W] = fastica(nM','numOfIC',2);
U = mean(nM,1);
tmp = bsxfun(@plus,A'*icasig,U');
close all
plot(pE)

%% get connected events
BLOB = ssdJ > .006;
CC = bwconncomp(BLOB);
%% filter front into wave front
wfront = front.*BLOB;
CCf = bwconncomp(wfront);
%% filter front2 
BLOB2 = ssd2J > .5*10^-3;
wfront2 = front2.*BLOB2;
%% filter front2 
BLOB0 = sJ > 2;
wfront0 = front0.*BLOB0;
%% for showing derivative
close all
tmpS = [];
for e = 1:size(dJ,3)
    tmpS(:,:,e) = imresize(sdJ(:,:,e),.5);
    e
end
%% show derivative
close all
for e = 1:size(sdJ,3)
    imshow(ssdJ(:,:,e),[])
    drawnow
end
%% explore cross section
close all
[r c v] = impixel(J(:,:,1),[]);
NP = norm([diff(r,1,1) diff(c,1,1)]);
r = linspace(r(1),r(2),round(NP));
c = linspace(c(1),c(2),round(NP));
WALL = zeros(300,round(NP));
sWALL = zeros(size(WALL));
for e = 1:300%size(J,3)
    imshow(J(:,:,e),[])
    hold on
    plot(r',c','g');
    WALL(e,:) = ba_interp2(sJ(:,:,e),r,c);
    sWALL(e,:) = ba_interp2(double(BLOB(:,:,e)),r,c);
    drawnow
    hold off
end
out = flattenMaskOverlay(bindVec(WALL),logical(sWALL),.3,'r');
figure;
imshow(out)
%% not needed
[d1 d2 d3] = gradient(sdJ);
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
    MSK = sdJ(:,:,e) > .01;
    imshow(sdJ(:,:,e).*MSK,[])
    %imshow(sdJ(:,:,e).*MSK,[])
    drawnow
end
%% show original
close all
for e = 1:3:size(sJ,3)
    imshow(sJ(:,:,e),[])
    title(num2str(e));
    drawnow
end
%% show movie with mask for wave fronts NEEDED to poke around for threshold
M3 = min(sJ(:));
M4 = max(sJ(:)-M3);
close all
for l = 1
    for e = 204%150:size(sJ,3)
        MSK = ssdJ(:,:,e) > .015;
        %MSK = ssdJ(:,:,e) > .008; % second leaf
        dB = bwboundaries(MSK);
        tmpF = logical(front(:,:,e));
        %tmpF = logical(wfront(:,:,e));
        out = flattenMaskOverlay((sJ(:,:,e)-M3)/M4,tmpF,.5,'b');
        tmpF = logical(front(:,:,e).*MSK);
        out = flattenMaskOverlay(out,tmpF,.5,'g');
        %out = cat(3,wv1(:,:,e),wv2(:,:,e),wv3(:,:,e));
        imshow(out,[]);
        %imshow((sJ(:,:,e)-M3)/M4,[])
        hold on
        for i = 1:numel(dB)
            %plot(dB{i}(:,2),dB{i}(:,1),'r')
        end
        hold off
        title(num2str(e));
        drawnow
        %MOV(:,:,:,e) = out;
    end
end
%%
for e = 1:size(MOV,4)
    sMOV(:,:,:,e) = imresize(MOV(:,:,:,e),.15);
end
%%
writerObj = VideoWriter('/home/nate/Downloads/leaf.avi');
open(writerObj);
close all
for e = 1:size(sMOV,4)
    imshow(MOV(:,:,:,e),[])
    drawnow
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);
%% plot all space time events
close all
figure
hold on
eSZ = [];
for k = 1:CCf.NumObjects
    [w1 w2 w3] = ind2sub(size(wfront),CCf.PixelIdxList{k});
    plot3(w1,w2,w3,'.')
    eSZ(k) = numel(w1);
end
%% measure speed
close all
%OLD THR = 70;
% Second movie is 150
[SP Ca_map TOT WVC decay stopMap conc DTM CMF] = mSpeed(CCf,size(wfront),70,ssdJ,sJ,wv1,wv2,wv3);

%% events over time
close all
TM = linspace(1,788,778);
for e = 1:size(DTM,1)
    f = normpdf(TM,DTM(e,2),DTM(e,3));
    f = bindVec(f)*DTM(e,4);
    vec = zeros(size(TM));
    vec(DTM(e,1)) = .2;
    plot(TM,f,'r')
    %plot(TM,vec,'b')
    hold on
    %waitforbuttonpress
end
%% plot fronts
close all
imshow(Ca_map,[]);
hold on
for e = 1:numel(WVC)
    plot(WVC{e}(:,2),WVC{e}(:,1),'r.')
end
%% mm
close all
[j1 j2] = ndgrid(linspace(-200,200,1000),linspace(-200,200,1000));
C = mvnpdf([j1(:) j2(:)],[0 0],[10 10]);
C = reshape(C,size(j1));
j1 = linspace(-200,200,1000);
C = normpdf(j1,0,5);
r = max(C);
for e = 1:100000
    C = C + .018*del2(C);
    plot(j1,C)
    axis([-200 200 0 r])
    drawnow
end
%% distribution
close all
xfac = 10/15;
[phat,pci] = mle(SP*xfac,'distribution','gam');
XS = linspace(0,20,100);
YS = pdf('Gamma',XS,phat(1),phat(2))
[yi xi] = ksdensity(SP*xfac);
yis = imfilter(yi,fspecial('average',[1 31]),'replicate');
plot(xi,yi,'b');
hold on
plot(xi,yis,'g');
plot(XS,YS,'r');
[J,fi] = max(YS);
title(num2str(XS(fi)));
figure;
plot(XS,YS,'r')
%% 
close all
didx = find(Ca_map~=0);
STACK = [];
D = [];
for t = 60:size(ssdJ,3)
    didx = find(wfront(:,:,t));
    LA = del2(sJ(:,:,t));
    dT = ssdJ(:,:,t);
    sX = LA(didx);
    sY = dT(didx);
    t
    %{
    plot(sX,sY,'.')
    %waitforbuttonpress
    axis equal
    title(num2str(t))
    corr(sX,sY)
    
    drawnow
    %}
    
     
    if ~isempty(didx)
        p = polyfit(sX,sY,1);
        D = [D;p(2)];
        if corr(sX,-sY) > .3
            STACK = [STACK;[sX sY]];
            plot(STACK(:,1),STACK(:,2),'.')
            drawnow
        end
    end
end
%% look at largest event
close all
[~,midx] = max(eSZ);
[JUNK,midx] = sort(eSZ);
midx = midx(end-10);
[w1 w2 w3] = ind2sub(size(wfront),CCf.PixelIdxList{midx});
plot3(w1,w2,w3,'.');
%% fill in area of influence
%% make the distance map for single events
singleBLOB = zeros(size(wfront));
singleBLOB(CCf.PixelIdxList{midx}) = 1;
OUTER = bwdist(singleBLOB);
INNER = bwdist(~singleBLOB);
DISTMAP = INNER - OUTER;
%% smooth distance map
sDISTMAP = smooth3(DISTMAP,'gaussian',[9 9 3]);
%% make surface ver0 of largest event
fv = isosurface(sDISTMAP,0);
[g1 g2 g3] = gradient(sDISTMAP);
%% get the normal to the curve - for tangent construction
n1 = zeros(size(sDISTMAP));
n2 = zeros(size(n1));
for e = 1:size(sDISTMAP,3)
    [n1(:,:,e) n2(:,:,e)] = gradient(sDISTMAP(:,:,e));
    e
end
%% show objects
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
            %{
            itype = 'linear';
            v1 = ba_interp2(n1(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),itype);
            v2 = ba_interp2(n2(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),itype);
            v3 = zeros(size(v1));
            %{
            V = (v1.^2 + v2.^2).^-.5;
            v1 = v1.*V;
            v2 = v2.*V;
            
            %}
            W1 = ba_interp2(g1(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),itype);
            W2 = ba_interp2(g2(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),itype);
            W3 = ba_interp2(g3(:,:,e),single(dB{k}(:,2)),single(dB{k}(:,1)),itype);
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
            vel = bsxfun(@times,C,C(:,3).^-1);
            
            quiver(dB{k}(:,2),dB{k}(:,1),vel(:,1),vel(:,2),'b')
            speed = sum(vel(:,1:2).*vel(:,1:2),2).^.5;
            
            cp = mean(dB{k},1);
           
            
            
            SP = [SP;speed];
            
            [MS loc] = max(speed);
            
            plot(dB{k}(loc,2),dB{k}(loc,1),'m*');
            
            quiver(dB{k}(:,2),dB{k}(:,1),v1,v2,'g');
            %text(cp(2),cp(1),num2str(max(speed),'% 10.2f'),'backgroundcolor','w');
            
            %quiver(dB{k}(:,2),dB{k}(:,1),t1,t2,'b');
            %}
        end
        
    end
    %figure(h2)
    %ksdensity(speed(:));
    %plot(w2(fidx),w1(fidx),'.')
    hold off
    drawnow
end
%% make the veins stick out
para.scales.value = 5;
para.resize.value = 1;
U = mean(sJ,3);
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
