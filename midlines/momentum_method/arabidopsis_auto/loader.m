%% load from model run
FilePath = '/mnt/scratch5/arabidopsis_model_13.09.09/';

FileList = {};
FileExt = {'mat'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% loop over saved data and load
parfor e = 1:numel(FileList)
    tic
    out(e) = loadData(FileList{e},20);
    toc
end
%% ingestigate the ME
for e = 1:numel(out)
    if ~isempty(out(e).ME)
        out(e).ME
        out(e).dataID
        e
    end
end
sele = 2601;
out(sele);
FilePath = '/mnt/scratch1/myLinks/arabidopsisGravitropism/baseLine/881/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
isolateRoots_overStack(SET{1},'/mnt/scratch5/arabidopsis_model_13.09.06/',1);
%%  loop over measurements SLOW
MK = [];
angle = [];
length = [];
growthRate = [];
frame = [];
comp = [];
for e = 1:numel(out)
    if isempty(out(e).ME) & (size(out(e).K,2) == 301)
        MK = cat(3,MK,out(e).K);
        angle = cat(2,angle,out(e).angle);
        length = cat(2,length,out(e).length);
        growthRate = cat(2,growthRate,out(e).growthRate);
        frame = cat(4,frame,out(e).frame);
        comp = cat(4,comp,permute(out(e).snipC,[3 2 1]));
    end
    fprintf(['done@' num2str(e) ':' num2str(numel(out)) '\n'])
end
%% loop over mesurements FAST
MK = zeros(200,301,numel(out));
angle = zeros(301,numel(out));
length = angle;
growthRate = angle;
frame = zeros(2,3,301,numel(out));
comp = zeros(2,301,300,numel(out));
for e = 1:numel(out)
    if isempty(out(e).ME) & (size(out(e).K,2) == 301)
        MK(:,:,e) = out(e).K;
        angle(:,e) = out(e).angle;
        length(:,e) = out(e).length;
        growthRate(:,e) = out(e).growthRate;
        frame(:,:,:,e) = out(e).frame;
        comp(:,:,:,e) = permute(out(e).snipC,[3 2 1]);
        rm(e) = 0;
    else
        rm(e) = 1;
    end
    fprintf(['done@' num2str(e) ':' num2str(numel(out)) '\n'])
end
%% remove those that did not load
ridx = find(rm);
angle(:,ridx) = [];
length(:,ridx) = [];
growthRate(:,ridx) = [];
frame(:,:,:,ridx) = [];
MK(:,:,ridx) = [];
comp(:,:,:,ridx) = [];
%out(ridx) = [];
%% clean
dA = diff(angle,1,1);
ridx = any(abs(dA) > 10*pi/180,1) | any(abs(growthRate) > 10,1) | any(abs(angle(1,:) - angle(end,:)) < 20*pi/180,1);
angle(:,ridx) = [];
length(:,ridx) = [];
growthRate(:,ridx) = [];
frame(:,:,:,ridx) = [];
MK(:,:,ridx) = [];
comp(:,:,:,ridx) = [];
%out(ridx) = [];
%% kmeans on angles
kidx = kmeans(angle',7);
UQ = unique(kidx);
f_ta = figure;
f_len = figure;
f_gr = figure;
for u = 1:numel(UQ)
    LEG{u} = num2str(UQ(u));
    
    figure;
    mesh(mean(MK(:,:,kidx==UQ(u)),3));
    view([0 90]);
    title(num2str(UQ(u)))
    
    figure(f_ta);
    hold all
    errorbar(mean(angle(:,kidx==UQ(u)),2),std(angle(:,kidx==UQ(u)),1,2)*sum(kidx==UQ(u))^-.5);
    
    figure(f_len);
    hold all
    errorbar(mean(length(:,kidx==UQ(u)),2),std(length(:,kidx==UQ(u)),1,2)*sum(kidx==UQ(u))^-.5);
    
    figure(f_gr);
    hold all
    errorbar(mean(growthRate(:,kidx==UQ(u)),2),std(growthRate(:,kidx==UQ(u)),1,2)*sum(kidx==UQ(u))^-.5);
end
figure(f_ta);
legend(LEG);
%% kmeans on growthRate
kidx = kmeans(mean(growthRate,1),100);
UQ = unique(kidx);
f_ta = figure;
f_len = figure;
f_gr = figure;
for u = 1:numel(UQ)
    LEG{u} = num2str(UQ(u));
    
    figure;
    mesh(mean(MK(:,:,kidx==UQ(u)),3));
    view([0 90]);
    title(num2str(UQ(u)))
    
    figure(f_ta);
    hold all
    errorbar(mean(angle(:,kidx==UQ(u)),2),std(angle(:,kidx==UQ(u)),1,2)*sum(kidx==UQ(u))^-.5);
    
    figure(f_len);
    hold all
    errorbar(mean(length(:,kidx==UQ(u)),2),std(length(:,kidx==UQ(u)),1,2)*sum(kidx==UQ(u))^-.5);
    
    figure(f_gr);
    hold all
    errorbar(mean(growthRate(:,kidx==UQ(u)),2),std(growthRate(:,kidx==UQ(u)),1,2)*sum(kidx==UQ(u))^-.5);
end
figure(f_ta);
legend(LEG);
%% model frame
norFrame = bsxfun(@minus,frame,frame(:,:,1,:));
T = myT(norFrame);
BV = T.decompose();
%% predict growthrate from curvuature
X = reshape(MK,[size(MK,1)*size(MK,2) size(MK,3)])';
Y = angle';
for L = 1:1

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,5);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);

        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);

        
        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}


        %{
        %%%%%%%%%%%%%%%
        % pls regression        
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);


end
%% view from above
for i = 1:size(k,1)
    plot(k(i,:),'r');
    hold on
    plot(Y(i,:),'b');
    hold off
    axis([0 301 0 pi/2]);
    drawnow
    pause(.5);
end
%% play with frames
delta = [[1 0 12];[0 1 13];[0 0 1]];
vec = [100 15];
vec = vec / norm(vec);
rvec = [-vec(2) vec(1)];
rot = [vec' rvec'];
%% ---------modeling frames BEGIN----------
%% create affine frames
for t = 1:size(frame,3)
    for e = 1:size(frame,4)
        nframe(:,:,t,e) = [frame(:,:,t,e);[0 0 1]];
    end
end
%% normalize to the first frame
for e = 1:size(nframe,4)
    initFrame(:,:,e) = squeeze(nframe(:,:,1,e));
    for t = 1:size(nframe,3)
        nframe(:,:,t,e) = inv(initFrame(:,:,e))*nframe(:,:,t,e);
        %nframe(:,:,t,e) = inv(mean(nframe(:,:,:,e),3))*nframe(:,:,t,e);
        %nframe(1:2,3,t,e) = nframe(1:2,3,t,e) - nframe(1:2,3,1,e);
    end
end
%% decompose
F = myT(nframe);
E = F.decompose();
%% ---------view FRAMES + COMP BEGIN----------
sidx = 700;
%displayModel({frame(:,:,:,sidx) frame(:,:,:,sidx) frame(:,:,:,sidx)},{comp(:,:,:,sidx) S.d(:,:,:,sidx) S2.d(:,:,:,sidx)});
%displayModel({frame(:,:,:,sidx) frame(:,:,:,sidx)},{comp(:,:,:,sidx) SMA(:,:,:,sidx)});
%displayModel({nframe(:,:,:,sidx) (:,:,:,sidx)},{comp(:,:,:,sidx) compS(:,:,:,sidx)});
displayModel({frame(:,:,:,sidx) frame(:,:,:,sidx)},{comp(:,:,:,sidx) compS(:,:,:,sidx)});
%% NOPE - - - - potenial function
F = squeeze(comp(:,:,:,sidx));
F = cat(1,F,ones([1,size(F,2),size(F,3)]));
dF = sum(F.*F,1).^-.5;
F = bsxfun(@times,F,dF);
F = permute(F,[2 3 1]);
mesh(cumsum(acos(sum(F(:,1:end-1,:).*F(:,2:end,:),3)),1))
%% ---------modeling comp BEGIN----------
X1 = myT(squeeze(comp(1,:,:,:)));
X2 = myT(squeeze(comp(2,:,:,:)));
E1 = X1.decompose();
E2 = X2.decompose();
%% select basis and test via sim
E1.setBasisNumber([4 4]);
E2.setBasisNumber([4 4]);
E1.generateBasisTensors();
E2.generateBasisTensors();
C1 = E1.project(X1);
C2 = E2.project(X2);
S1 = E1.construct(C1);
S2 = E2.construct(C2);
SMA = cat(1,shiftdim(S1.d,-1),shiftdim(S2.d,-1));
%% view 
plot(X2.d(200,:,10));
hold all;
plot(S2.d(200,:,10));
%% group X1 and X2 together
X = myT(squeeze(comp(:,:,:,:)));
E = X.decompose();
%% select on X1 and X2 together
E.setBasisNumber([2 4 4]);
E.generateBasisTensors();
C = E.project(X);
S = E.construct(C);
%% view 
plot(squeeze(X.d(2,301,:,100)));
hold all;
plot(squeeze(S.d(2,301,:,100)));
%% decompose further
data = C.d;
[dS dC dU dE dL dERR dLAM] = PCA_FIT_FULL(data,10);
dS = myT(dS);
S2 = E.construct(dS);
%% plot ov kmeans
figure;
CM = [];
uG = [];
kidx = kmeans(mean(growthRate,1),20);
kidx = kmeans(dC,20);
UQ = unique(kidx);

for u = 1:numel(UQ)
    fidx = find(kidx==UQ(u));
    uG(u) = mean(mean(growthRate(:,fidx,1)));
    CM(u,:) = mean(dC(kidx==UQ(u),:));
end
[J sidx] = sort(uG);
CM = CM(sidx,:);
%%
for u = 1:numel(UQ)
    plot3(dC(kidx==UQ(u),1),dC(kidx==UQ(u),2),dC(kidx==UQ(u),3),'.');
    hold all
    LEG{u} = num2str(u);
end
%%
plot3(CM(:,1),CM(:,2),CM(:,3),'k');
legend(LEG)
%%
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL(growthRate',4);
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL(angle',4);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try bundle model direct sum of X1 X2;
X1 = squeeze(comp(1,:,:,:));
X2 = squeeze(comp(2,:,:,:));
X1 = permute(X1,[2 1 3]);
X2 = permute(X2,[2 1 3]);
sz1 = size(X1);
sz2 = size(X2);
X1 = reshape(X1,[size(X1,1) size(X1,2)*size(X1,3)]);
X2 = reshape(X2,[size(X2,1) size(X2,2)*size(X2,3)]);
[X1S X1C X1U X1E X1L X1ERR X1LAM] = PCA_FIT_FULL(X1',4);
[X2S X2C X2U X2E X2L X2ERR X2LAM] = PCA_FIT_FULL(X2',4);
rX1C = reshape(X1C',[4 sz1(2) sz1(3)]);
rX2C = reshape(X2C',[4 sz2(2) sz2(3)]);
X1 = reshape(X1,sz1);
X1S = reshape(X1S',sz1);
X2 = reshape(X2,sz2);
X2S = reshape(X2S',sz2);
X1 = ipermute(X1,[2 1 3]);
X2 = ipermute(X2,[2 1 3]);
X1S = ipermute(X1S,[2 1 3]);
X2S = ipermute(X2S,[2 1 3]);
compS = cat(1,shiftdim(X1S,-1),shiftdim(X2S,-1));
testS = cat(1,shiftdim(X1,-1),shiftdim(X2,-1));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try pca on finite K
T = diff(comp,1,3);
lT = sum(T.*T,1).^.5;
T = bsxfun(@times,T,lT.^-1);
angle = squeeze(atan2(T(2,:,:,:),T(1,:,:,:)));
K = diff(angle,1,2);
%quiver(squeeze(comp(1,200,1:end-1,1)),squeeze(comp(2,200,1:end-1,1)),squeeze(T(1,200,:,1)),squeeze(T(2,200,:,1)))
X1 = angle;

X1 = permute(X1,[2 1 3]);
sz1 = size(X1);
X1 = reshape(X1,[size(X1,1) size(X1,2)*size(X1,3)]);
[XS XC XU XE XL XERR XLAM] = PCA_FIT_FULL(X1',3);

% restore orginal data
XS = reshape(XS',sz1);


XC = reshape(XC',[3 sz1(2:end)]);
XC = permute(XC,[2 1 3]);
%% plot direct sum bundle model of K
for e = 1:10%size(rXC,3)   
    plot3(XC(1:end,1,e),XC(1:end,2,e),XC(1:end,3,e));
    hold on
    plot3(XC(1,1,e),XC(1,2,e),XC(1,3,e),'r*');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bundle direct sum X1,X2, affine frame
X1 = squeeze(comp(1,:,:,:));
X2 = squeeze(comp(2,:,:,:));
fsz = size(nframe);
X3 = reshape(nframe,[fsz(1)*fsz(2) fsz(3) fsz(4)]);
X1 = permute(X1,[2 1 3]);
X2 = permute(X2,[2 1 3]);
sz1 = size(X1);
sz2 = size(X2);
sz3 = size(X3);

X1 = reshape(X1,[size(X1,1) size(X1,2)*size(X1,3)]);
X2 = reshape(X2,[size(X2,1) size(X2,2)*size(X2,3)]);
X3 = reshape(X3,[size(X3,1) size(X3,2)*size(X3,3)]);

D = [X1' X2' X3'];
sD = std(D,1,1);
uD = mean(D,1);
D = bsxfun(@minus,D,uD);
sD(sD==0) = 1;
D = bsxfun(@times,D,sD.^-1);
[XS XC XU XE XL XERR XLAM] = PCA_FIT_FULL(D,6);
XS = bsxfun(@times,XS,sD);
XS = bsxfun(@plus,XS,uD);
% restore orginal data
X1 = reshape(X1,sz1);
X2 = reshape(X2,sz2);
X3 = reshape(X3,sz3);

X1S = XS(:,1:sz1(1));
X2S = XS(:,sz1(1)+1:sz1(1)+sz2(1));
X3S = XS(:,sz1(1)+sz2(1)+1:end);
X1S = reshape(X1S',sz1);
X2S = reshape(X2S',sz2);
X3S = reshape(X3S',sz3);
X3S = reshape(X3S,fsz);
X1 = ipermute(X1,[2 1 3]);
X2 = ipermute(X2,[2 1 3]);
X1S = ipermute(X1S,[2 1 3]);
X2S = ipermute(X2S,[2 1 3]);
compS = cat(1,shiftdim(X1S,-1),shiftdim(X2S,-1));
frameS = X3S;
%%
rXC = reshape(XC',[6 sz1(2) sz1(3)]);
rXC = ipermute(rXC,[2 1 3]);
%% plot direct sum bundle model
for e = 1:100%size(rXC,3)   
    plot3(rXC(1:end,1,e),rXC(1:end,2,e),rXC(1:end,3,e));
    hold on
    plot3(rXC(1,1,e),rXC(1,2,e),rXC(1,3,e),'r*');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bundle model direct sum over time HERE IS BEST ONE SO FAR
X1 = squeeze(comp(1,:,:,:));
X2 = squeeze(comp(2,:,:,:));
X1 = permute(X1,[2 1 3]);
X2 = permute(X2,[2 1 3]);
sz1 = size(X1);
sz2 = size(X2);

X1 = reshape(X1,[size(X1,1) size(X1,2)*size(X1,3)]);
X2 = reshape(X2,[size(X2,1) size(X2,2)*size(X2,3)]);

[XS XC XU XE XL XERR XLAM] = PCA_FIT_FULL([X1' X2'],3);

X1S = XS(:,1:sz1(1));
X2S = XS(:,sz1(1)+1:end);

rXC = reshape(XC',[3 sz1(2) sz1(3)]);
rXC = ipermute(rXC,[2 1 3]);


X1 = reshape(X1,sz1);
X1S = reshape(X1S',sz1);
X2 = reshape(X2,sz2);
X2S = reshape(X2S',sz2);
X1 = ipermute(X1,[2 1 3]);
X2 = ipermute(X2,[2 1 3]);
X1S = ipermute(X1S,[2 1 3]);
X2S = ipermute(X2S,[2 1 3]);
compS = cat(1,shiftdim(X1S,-1),shiftdim(X2S,-1));
%% plot direct sum bundle model
for e = 1:5%size(rXC,3)   
    plot3(rXC(:,1,e),rXC(:,2,e),rXC(:,3,e));
    hold on
    plot3(rXC(1,1,e),rXC(1,2,e),rXC(1,3,e),'r*');
end
%% pca for denoise 
SMOOTH_DIM = 7;
szC = size(rXC);
rXC = reshape(rXC,[szC(1) szC(2)*szC(3)]);
[SS SC SU SE SL SERR SLAM] = PCA_FIT_FULL(rXC',SMOOTH_DIM);
SS = reshape(SS',szC);
% for testing SMOOTH_DIM
test = SS;
test = permute(test,[2 1 3]);
test = reshape(test,[size(test,1) size(test,2)*size(test,3)]);
test_mani = PCA_BKPROJ(test',XE,XU);

test_X1S = test_mani(:,1:sz1(1));
test_X2S = test_mani(:,sz1(1)+1:end);
test_X1S = reshape(test_X1S',sz1);
test_X2S = reshape(test_X2S',sz2);
test_X1S = ipermute(test_X1S,[2 1 3]);
test_X2S = ipermute(test_X2S,[2 1 3]);
test_comp = cat(1,shiftdim(test_X1S,-1),shiftdim(test_X2S,-1));


SC = reshape(SC',[SMOOTH_DIM szC(2:end)]);
rXC = reshape(rXC,szC);
%% last decompose of core tensors 7 3
z = size(SC);
SC = reshape(SC,[z(1)*z(2) z(3)]);
[fS fC fU fE fL fERR fLAM] = PCA_FIT_FULL(SC',5);
SC = reshape(SC,[z(1) z(2) z(3)]);
%% plot direct sum bundle model and denoise
for e = 1:100%size(rXC,3)   
    %plot3(rXC(:,1,e),rXC(:,2,e),rXC(:,3,e),'k');
    hold on
    %plot3(rXC(1,1,e),rXC(1,2,e),rXC(1,3,e),'r*');
    %plot3(rXC(end,1,e),rXC(end,2,e),rXC(end,3,e),'b*');
    
    plot3(SS(:,1,e),SS(:,2,e),SS(:,3,e),'k');
    plot3(SS(1,1,e),SS(1,2,e),SS(1,3,e),'r*');
    plot3(SS(end,1,e),SS(end,2,e),SS(end,3,e),'b*');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bundle model direct sum WHOLE
X1 = squeeze(comp(1,:,:,:));
X2 = squeeze(comp(2,:,:,:));
X1 = permute(X1,[2 1 3]);
X2 = permute(X2,[2 1 3]);
sz1 = size(X1);
sz2 = size(X2);

X1 = reshape(X1,[size(X1,1)*size(X1,2) size(X1,3)]);
X2 = reshape(X2,[size(X2,1)*size(X2,2) size(X2,3)]);

[XS XC XU XE XL XERR XLAM] = PCA_FIT_FULL([X1' X2'],3);

X1S = XS(:,1:sz1(1));
X2S = XS(:,sz1(1)+1:end);
rXC = reshape(XC',[3 sz1(2) sz1(3)]);
%% BEGIN **** projective kmeans modeling
% choose to decompose comp together
X = myT(squeeze(comp(:,:,:,:)));
E = X.decompose();
E.setBasisNumber([2 4 4]);
E.generateBasisTensors();
C = E.project(X);
S = E.construct(C);
% project core tensor C into 3
[coreS coreC coreU coreE coreL coreERR coreLAM] = PCA_FIT_FULL(C.d,5);
coreS = myT(coreS);
S = E.construct(C);
%%

figure;
CM = [];
uG = [];
kidx = kmeans(mean(growthRate,1),20);
UQ = unique(kidx);
for u = 1:numel(UQ)
    fidx = find(kidx==UQ(u));
    uG(u) = mean(mean(growthRate(:,fidx,1)));
    CM(u,:) = mean(coreC(kidx==UQ(u),:));
end
[J sidx] = sort(uG);
CM = CM(sidx,:);
plot3(coreC(:,1),coreC(:,2),coreC(:,3),'.','MarkerSize',1);
hold on;
plot3(CM(:,1),CM(:,2),CM(:,3),'k')
%%
dC = diff(CM,1,1);
dC = sum(dC.*dC,2).^.5;
dC = cumsum([0;dC]);
CM = interp1(dC,CM,linspace(0,dC(end),100));
pkmeans(CM,coreC,10);
%% BEGIN AGAIN
T = diff(comp,1,3);
lT = sum(T.*T,1).^.5;
T = bsxfun(@times,T,lT.^-1);
angle = squeeze(atan2(T(2,:,:,:),T(1,:,:,:)));
%{
for e = 1:size(angle,3)
    mesh(angle(:,:,e))
    view([0 90]);
    drawnow
pause(.3)
end
%}
K = diff(angle,1,2);
X1 = K;
X1 = permute(X1,[2 1 3]);
sz1 = size(X1);
X1 = reshape(X1,[size(X1,1) size(X1,2)*size(X1,3)]);
[XS XC XU XE XL XERR XLAM] = PCA_FIT_FULL(X1',2);
K = reshape(XS',sz1);
[P,A,W] = fastica(XS');


dE = sum(K,1);
F2 = (moment(K,2,1));
FEA = cat(1,dE,F2);
FEA = reshape(FEA,[size(FEA,1) size(FEA,2)*size(FEA,3)]);


dE_comp = permute(comp,[1 3 2 4]);
szE = size(dE_comp);
dE_comp = reshape(dE_comp,[szE(1) szE(2) szE(3)*szE(4)]);
[J sidx] = sortrows(FEA');
dE_comp = dE_comp(:,:,sidx);
%%
for e = 1:300:size(dE_comp,3)
    %plot(dE_comp(1,:,e),dE_comp(2,:,e));
    plot(dE_comp(2,:,e))
    axis equal
    drawnow
end
%% AND AGAIN
X1 = squeeze(comp(1,:,:,:));
X2 = squeeze(comp(2,:,:,:));

sz1 = size(X1);
sz2 = size(X2);
sznf = size(nframe);
sz = sz1;
X1 = reshape(X1,[size(X1,1)*size(X1,2) size(X1,3)]);
X2 = reshape(X2,[size(X2,1)*size(X2,2) size(X2,3)]);
X3 = reshape(nframe,[sznf(1)*sznf(2)*sznf(3) sznf(4)]);

[XS XC XU XE XL XERR XLAM] = PCA_FIT_FULL([X1' X2'],3);

X1  = reshape(X1,sz1);
X2  = reshape(X2,sz2);

X1S = XS(:,1:sz1(1)*sz1(2))';
X2S = XS(:,(sz1(1)*sz1(2)+1):end)';

X1S  = reshape(X1S,sz1);
X2S  = reshape(X2S,sz2);

compS = cat(1,shiftdim(X1S,-1),shiftdim(X2S,-1));
%% view and again
figure
plot3(XC(:,1),XC(:,2),XC(:,3),'.','MarkerSize',5);
hold on;
%% kmeans n growthrate on and again
CM = [];
uG = [];
kidx = kmeans(mean(growthRate,1),12);
UQ = unique(kidx);

for u = 1:numel(UQ)
    fidx = find(kidx==UQ(u));
    uG(u) = mean(mean(growthRate(:,fidx,1)));
    CM(u,:) = mean(XC(kidx==UQ(u),:));
end
[J sidx] = sort(uG);
CM = CM(sidx,:);
plot3(CM(:,1),CM(:,2),CM(:,3),'k')
dL = diff(CM,1,1);
dL = sum(dL.*dL,2).^.5;
L = cumsum([0;dL]);
fun = interp1(L,CM,linspace(0,L(end),1000),'spline');
plot3(fun(:,1),fun(:,2),fun(:,3),'b');
[U E lam CS] = pkmeans(fun,XC,51);
lam= lam.^.5;
%% align vectors
for e = 2:size(E,3)
    if sign(E(:,1,e)'*E(:,1,e-1)) < 0
        E(:,1,e) = -E(:,1,e);
    end
    
    if sign(E(:,2,e)'*E(:,2,e-1)) < 0
        E(:,2,e) = -E(:,2,e);
    end
end
%% plot unwrapped space
figure;
plot3(CS(:,1),CS(:,2),CS(:,3),'.')
%% plot model space
model = figure;
plot3(XC(:,1),XC(:,2),XC(:,3),'.','MarkerSize',5);
hold on;
mag = 500;
for e = 1:10:size(U,1)
    quiver3(U(e,1),U(e,2),U(e,3),E(1,1,e),E(2,1,e),E(3,1,e),lam(1),'r');
    %quiver3(U(e,1),U(e,2),U(e,3),-E(1,1,e),-E(2,1,e),-E(3,1,e),mag,'r');
    quiver3(U(e,1),U(e,2),U(e,3),E(1,2,e),E(2,2,e),E(3,2,e),lam(2),'g');
    %quiver3(U(e,1),U(e,2),U(e,3),-E(1,2,e),-E(2,2,e),-E(3,2,e),mag,'g');
end
%C1 = U + 500*permute(squeeze(E(:,1,:)),[2 1 3]);
%plot3(C1(:,1),C1(:,2),C1(:,3),'k')
%% measure model along parameters
c = C1(1,:);
% measure along center curve
hk = figure;
for e = 1:10:size(U,1)
    K = measureModel(XU,XE,U(e,:),sz);
    figure(hk);
    mesh(K)
    view([0 90]);
    drawnow
end
%% measure along E1 @ some pos
magy = 2;
pos = size(U,1)/2;
pos = size(U,1)/4;
pos = 3*size(U,1)/4;
pos = size(U,1);
pos = round(pos);
svec = 1;
p1 = linspace(-magy*lam(svec,pos),magy*lam(svec,pos),30);
for e = 1:numel(p1)
    vecT(e,:) = U(pos,:)+p1(e)*E(:,svec,pos)';
    K = measureModel(XU,XE,vecT(e,:),sz);
    figure(hk);
    mesh(-K)
    view([0 90]);
    drawnow
end
figure(model);
plot3(vecT(:,1),vecT(:,2),vecT(:,3),'m','LineWidth',3);
%%
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL(growthRate',4);
[aS aC aU aE aL aERR aLAM] = PCA_FIT_FULL(angle',4);
modelViewer(XU,XE,sz,U,E,lam,XC,gS',aS',frame);

%% frame modeling
for tr = 1:size(frame,4)
    loc_stack = squeeze(frame(:,3,:,tr))';
    [fS fC fU(tr,:) fE(:,:,tr) fL fERR fLAM] = PCA_FIT_FULL(loc_stack,2);
    dv = mean(diff(loc_stack,1,1));
    if dv*fE(:,1,tr) < 0
        fE(:,1,tr) = -fE(:,1,tr);
    end
    fE(1,2,tr) = -fE(2,1,tr);
    fE(2,2,tr) = fE(1,1,tr);
    fC = PCA_REPROJ(loc_stack,fE(:,:,tr),fU(tr,:));
    nframe(:,:,:,tr) = frame(:,:,:,tr);
    nframe(:,3,:,tr) = fC';
end
%%
sznf = size(nframe);
rnframe = reshape(nframe,[sznf(1)*sznf(2)*sznf(3) sznf(4)]);
[nfS nfC nfU nfE nfL nfERR nfLAM] = PCA_FIT_FULL(rnframe',3);
fS = reshape(fS',sznf);
%%

for i = 1:size(nframe,4)
    for t = 2:size(nframe,3)
        L(t,i) = norm(nframe(1:2,3,t,i));
        nframe(1:2,3,t,i) = nframe(1:2,3,t,i)/L(t,i);
    end
end
%%
sznf = size(nframe);
rnframe = reshape(nframe,[sznf(1)*sznf(2)*sznf(3) sznf(4)]);
[fS fC fU fE fL fERR fLAM] = PCA_FIT_FULL(rnframe',6);
fS = reshape(fS',sznf);
%% view frame scatter
figure;
plot3(nfC(:,1),nfC(:,2),nfC(:,3),'.')
%%
for i = 1:size(nframe,4)
    for t = 2:size(nframe,3)
        fS(1:2,3,t,i) = fS(1:2,3,t,i)*L(t,i);
        nframe(1:2,3,t,i) = nframe(1:2,3,t,i)*L(t,i);
    end
end
%%
figure;
tr = 200;
for t = 1:10:size(fS,3)
    quiver3(fS(1,3,t,tr),fS(2,3,t,tr),t,fS(1,2,t,tr),fS(2,2,t,tr),0,'r');
    hold on
    quiver3(fS(1,3,t,tr),fS(2,3,t,tr),t,fS(1,1,t,tr),fS(2,1,t,tr),0,'r');

    quiver3(nframe(1,3,t,tr),nframe(2,3,t,tr),t,nframe(1,2,t,tr),nframe(2,2,t,tr),0,'k');
    quiver3(nframe(1,3,t,tr),nframe(2,3,t,tr),t,nframe(1,1,t,tr),nframe(2,1,t,tr),0,'k');
end





