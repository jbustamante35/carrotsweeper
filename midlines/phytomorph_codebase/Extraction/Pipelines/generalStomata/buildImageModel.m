function [imageGradeFunc] = buildImageModel(imageData,domainMask,osz,modelDIM)

%% BUILD MODELS
close all


gidxI = find(domainMask==1);
gidxE = find(domainMask==0);

ecach = [];
icach = [];


%d = cat(3,imageData,flipdim(imageData,1),flipdim(imageData,2),flipdim(flipdim(imageData,1),2));

for e = 1:size(imageData,1)
    tmp = reshape(imageData(e,:),osz);
    %tmp = circshift(tmp,100,2);
    tmp = flipdim(tmp,2);
    d = [imageData(e,:);tmp(:)'];

    [d,do] = patchFlipOp(reshape(imageData(e,:),osz));
    d = d';
    m1 = bsxfun(@times,do,domainMask);
    idx = nchoosek(1:4,2);
    for m = 1:size(idx,1)
        dV = m1(:,:,idx(m,1)) - m1(:,:,idx(m,2));
        Tmpmea1(m) = norm(dV(:));
    end
    mea1(e) = mean(Tmpmea1);
    %d = d';
    d = imageData(e,:);
    d = zscore(d,1,2);
    icach = [icach;d(:,gidxI)];
    ecach = [ecach;d(:,gidxE)];
%{
    icach(e,:) = imageData(e,gidxI);
    ecach(e,:) = imageData(e,gidxE);
    %}
    

    %icach2(e,:) = cach2(e,gidxI);
    %ecach2(e,:) = cach2(e,gidxE);
end
%icach = bsxfun(@minus,icach,mean(icach,2));
%icach = zscore(icach,1,2);
[icU,icE,icL]= PCA_FIT_FULLws(icach,modelDIM(1));
%{
nicU = icU/norm(icU);
nicE = icE - nicU'*(nicU*icE);
for e = 1:size(nicE,2)
    nicE(:,e) = nicE(:,e)/norm(nicE(:,e));
end
icE = nicE;
%}
%ecach = zscore(ecach,1,2);
[ecU,ecE,ecL]= PCA_FIT_FULLws(ecach,modelDIM(2));
%{
necU = ecU/norm(ecU);
necE = ecE - necU'*(necU*ecE);
for e = 1:size(necE,2)
    necE(:,e) = necE(:,e)/norm(necE(:,e));
end
ecE = necE;
%}

[icC]= PCA_REPROJ(icach,icE,icU);
[ecC]= PCA_REPROJ(ecach,ecE,ecU);
figure;
[iSim]= PCA_BKPROJ(icC,icE,icU);
[eSim]= PCA_BKPROJ(ecC,ecE,ecU);
icCS = icC*diag(diag(icL).^-1);
ecCS = ecC*diag(diag(ecL).^-1);



idelta = (iSim - icach);
[idU,idE,idL]= PCA_FIT_FULLws(idelta,2);
[idC]= PCA_REPROJ(idelta,idE,idU);
ksId = [idC];
% ksdensity plots
for e = 1:size(ksId,2)
    figure
    [ipdffD(e,:),ixiD(e,:)] = ksdensity(ksId(:,e),'Bandwidth',.25*std(ksId(:,e)));
    ipdffD(e,:) = ipdffD(e,:)/sum(ipdffD(e,:));
    plot(ixiD(e,:),ipdffD(e,:));
    title(['PC:' num2str(e) '<--internal']);
end

%{
for e1 = 1:size(icC,1)
    for e2 = 1:size(icC,1)
        iDIST(e1,e2) = norm(icC(e1,:) - icC(e2,:));
    end
end

figure;
path = 3*[cos(linspace(-pi,pi,40))' sin(linspace(-pi,pi,40))'];
[piSim]= PCA_BKPROJ([0 0],icE,icU);
Z = [];
for p = 1:size(path,1)
   [peSim]= PCA_BKPROJ(path(p,:),ecE,ecU); 
    z = zeros(size(Gmask));
    %z(gidxI) = piSim;
    z(gidxE) = ecE(:,2);
    imshow(z,[-1 1]);
    dimageDatanow
    Z = [Z z];
    %waitforbuttonpress
end

figure;
imshow(Z,[]);
%}
if size(icCS,2) > 1
    figure;
    plot(icCS(:,1),icCS(:,2),'.')
    figure;
    plot(ecCS(:,1),ecCS(:,2),'.');
end
%hold on
%plot(path(:,1),path(:,2),'r')

%waitforbuttonpress

figure
iErr = sum((iSim - icach).^2,2).^.5;
eErr = sum((eSim - ecach).^2,2).^.5;
iPDF_u = [zeros(1,size(icE,2)) mean(iErr)];
iPDF_c = diag([diag(icL);std(iErr)^2]);
ePDF_u = [zeros(1,size(ecE,2)) mean(eErr)];
ePDF_c = diag([diag(ecL);std(eErr)^2]);
maskI = zeros(osz);
maskI(gidxI) = icU;
imshow(maskI,[]);
maskE = zeros(osz);
maskE(gidxE) = ecU;
imshow(maskE,[]);
maskT = zeros(size(maskI));
maskT = maskI + maskE;
imshow(maskT,[]);
maskW = reshape(mean(imageData,1),osz);
imshow(maskW,[]);
MM = [];
for d = 1:size(icE,2)
    tM = [];
    sweep = linspace(-icL(d,d),icL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(icE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,icE,icU);
        tmpM = zeros(size(maskI));
        tmpM(gidxI) = tmpD;
        tM = [tM tmpM];
    end
    MM = [MM;tM];
end
imshow(MM,[]);
figure;
MM = [];
for d = 1:size(ecE,2)
    tM = [];
    sweep = linspace(-ecL(d,d),ecL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(ecE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,ecE,ecU);
        tmpM = zeros(size(maskI));
        tmpM(gidxE) = tmpD;
        tM = [tM tmpM];
    end
    MM = [MM;tM];
end
imshow(MM,[]);
% regress
rIE = icC\ecC;
rEI = ecC\icC;
ksI = [icC iErr];
ksE = [ecC eErr];
%
clear ipdff epdff
% ksdensity plots
for e = 1:size(ksI,2)
    figure
    [ipdff(e,:),ixi(e,:)] = ksdensity(ksI(:,e),'Bandwidth',.25*std(ksI(:,e)));
    ipdff(e,:) = ipdff(e,:)/sum(ipdff(e,:));
    plot(ixi(e,:),ipdff(e,:));
    title(['PC:' num2str(e) '<--internal']);
end
for e = 1:size(ksE,2)
    figure
    [epdff(e,:),exi(e,:)] = ksdensity(ksE(:,e),'Bandwidth',.5*std(ksE(:,e)));
    epdff(e,:) = epdff(e,:)/sum(epdff(e,:));
    plot(exi(e,:),epdff(e,:));
    title(['PC:' num2str(e) '<--external']);
end
clear Xi2 Xe2

nPTS = 30;
Mxi = [(1:nPTS)' (1:nPTS)'];
Myi = [linspace(min(ksI(:,1)),max(ksI(:,1)),nPTS)' linspace(min(ksI(:,2)),max(ksI(:,2)),nPTS)'];
if size(ksI,2) == 2
    Myi = [linspace(min(ksI(:,1)),max(ksI(:,1)),nPTS)'];
end
MXui = mean(Mxi,1);
MYui = mean(Myi,1);
Mxi = bsxfun(@minus,Mxi,MXui);
Myi = bsxfun(@minus,Myi,MYui);
MAP1(1,1) = Myi(:,1)\Mxi(:,1);
if size(ksI,2) > 2
    MAP1(2,2) = Myi(:,2)\Mxi(:,2);
else
    MAP1(2,2) = 0;
end
MAP1(1,2) = 0;
MAP1(2,1) = 0;
MAP1(1,3) = MXui(1);
MAP1(2,3) = MXui(2);
MAP1(3,3) = 1;
if numel(MYui) ~= 1
    MAPYY1 = [eye(2) -MYui'];
else
    MAPYY1 = [eye(2) -[MYui 0]'];
end
MAPYY1(3,3) = 1;
MAP1 = MAP1*MAPYY1;
MAP1(3,:) = [];
Myi = bsxfun(@plus,Myi,MYui);
if size(Myi,2) == 1
    test1 = (MAP1*[Myi zeros(size(Myi)) ones(size(Myi,1),1)]')';
else
    test1 = (MAP1*[Myi ones(size(Myi,1),1)]')';
end


tp = [Myi(1,1) 0];
testTP = (MAP1*[tp ones(size(tp,1),1)]')';
figure;
plot(test1(:,1),Mxi(:,1),'.')
figure;
plot(test1(:,2),Mxi(:,2),'.')

Mxe = [(1:nPTS)' (1:nPTS)'];
Mye = [linspace(min(ksE(:,1)),max(ksE(:,1)),nPTS)' linspace(min(ksE(:,2)),max(ksE(:,2)),nPTS)'];
MXue = mean(Mxe,1);
MYue = mean(Mye,1);
Mxe = bsxfun(@minus,Mxe,MXue);
Mye = bsxfun(@minus,Mye,MYue);
MAP2(1,1) = Mye(:,1)\Mxe(:,1);
MAP2(2,2) = Mye(:,2)\Mxe(:,2);
MAP2(1,2) = 0;
MAP2(2,1) = 0;
MAP2(1,3) = MXue(1);
MAP2(2,3) = MXue(2);
MAP2(3,3) = 1;
MAPYY1 = [eye(2) -MYue'];
MAPYY1(3,3) = 1;
MAP2 = MAP2*MAPYY1;
MAP2(3,:) = [];
Mye = bsxfun(@plus,Mye,MYue);
test2 = (MAP2*[Mye ones(size(Mye,1),1)]')';
figure;
plot(test2(:,1),Mxe(:,1),'.')


[Xi2(:,:,1),Xi2(:,:,2)] = ndgrid(linspace(min(ksI(:,1)),max(ksI(:,1)),nPTS),linspace(min(ksI(:,2)),max(ksI(:,2)),nPTS));
XiX = Xi2(:,:,1);
XiY = Xi2(:,:,2);
[Ki2,Xi2] = ksdensity(ksI(:,1:2),[XiX(:) XiY(:)],'Bandwidth',.5*std(ksI(:,1:2),1,1));
Xi2 = reshape(Xi2,[nPTS nPTS 2]);
Ki2 = reshape(Ki2,[nPTS nPTS]);
[Xe2(:,:,1),Xe2(:,:,2)] = ndgrid(linspace(min(ksE(:,1)),max(ksE(:,1)),nPTS),linspace(min(ksE(:,2)),max(ksE(:,2)),nPTS));
XeX = Xe2(:,:,1);
XeY = Xe2(:,:,2);
[Ke2,Xe2] = ksdensity(ksE(:,1:2),[XeX(:) XeY(:)],'Bandwidth',.5*std(ksE(:,1:2),1,1));
Xe2 = reshape(Xe2,[nPTS nPTS 2]);
Ke2 = reshape(Ke2,[nPTS nPTS]);
figure;
mesh(Xi2(:,:,1),Xi2(:,:,2),Ki2);
figure;
mesh(Xe2(:,:,1),Xe2(:,:,2),Ke2);
figure;
plot(sum(Ki2,2))
figure;
plot(ipdff(1,:))
%
% store PDF data
clear myiPDF myePDF 
myiPDF.X = ixi;
myiPDF.F = ipdff;
myiPDF.XD = ixiD;
myiPDF.FD = ipdffD;
myePDF.X = exi;
myePDF.F = epdff;
myiPDF.X2Xi = Xi2; 
myiPDF.X2Yi = Ki2/sum(Ki2(:));
myePDF.X2Xe = Xe2; 
myePDF.X2Ye = Ke2/sum(Ke2(:));
myiPDF.MAP1 = MAP1;
myePDF.MAP2 = MAP2;


[mf1,mfx] = ksdensity(mea1);
mf1 = mf1 / sum(mf1);
moM = [mfx;mf1];

% regression plots
toScatterE = (icC*rIE);
for r = 1:size(toScatterE,2)
    figure;
    plot(ecC(:,r),toScatterE(:,r),'k.')
    hold on
    plot(linspace(min(ecC(:,r)),max(ecC(:,r)),10),linspace(min(ecC(:,r)),max(ecC(:,r)),10),'r');
end
% regression plots
toScatterI = (ecC*rEI);
for r = 1:size(toScatterI,2)
    figure;
    plot(icC(:,r),toScatterI(:,r),'k.')
    hold on
    plot(linspace(min(icC(:,r)),max(icC(:,r)),10),linspace(min(icC(:,r)),max(icC(:,r)),10),'r');
end
figure;
MM = [];
for d = 1:size(icE,2)
    tM = [];
    sweep = linspace(-icL(d,d),icL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(icE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,icE,icU);
        tmpM = zeros(size(maskI));
        tmpM(gidxI) = tmpD;
    
        tmpE = (tmpC*rIE);
        tmpE = PCA_BKPROJ(tmpE,ecE,ecU);
        tmpME = zeros(size(maskI));
        tmpME(gidxE) = tmpE;
        tmpM = tmpM + tmpME;

        tM = [tM tmpM];
    end
    MM = [MM;tM];
end

imshow(MM,[]);
title('internal->external');
drawnow
figure;
MM = [];
for d = 1:size(ecE,2)
    tM = [];
    sweep = linspace(-ecL(d,d),ecL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(ecE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,ecE,ecU);
        tmpM = zeros(size(maskI));
        tmpM(gidxE) = tmpD;
    
        tmpI = (tmpC*rEI);
        tmpI = PCA_BKPROJ(tmpI,icE,icU);
        tmpMI = zeros(size(maskI));
        tmpMI(gidxI) = tmpI;
        tmpM = tmpM + tmpMI;

        tM = [tM tmpM];
    end
    MM = [MM;tM];
end

imshow(MM,[]);
title('external->internal');
% pack into model function
imageGradeFunc = @(X)issueUnitGrades(X,gidxI,icE,icU,iPDF_u,iPDF_c,gidxE,ecE,ecU,ePDF_u,ePDF_c,rEI,rIE,idU,idE,myiPDF,myePDF,moM,domainMask);

end