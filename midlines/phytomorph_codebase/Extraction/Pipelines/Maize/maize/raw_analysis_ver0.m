%% 
oGravi = bsxfun(@minus,mGravi,mGravi(:,1));
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL(mGravi,4);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(mSpec,500);
close all
vals = [];
opt = statset('UseParallel','always');
X = sC;
Y = gC;
ITR = 5;
KF = 5;
for tr = 1:ITR
    for n = 1:20
        func = @(X1,Y1,X2,Y2)myCV_ver0(X1,Y1,X2,Y2,n);
        func = @(TR,TE)myCV_ver1_maize(TR,TE,n);
        vals(:,:,n,tr) = crossval(func,[X,Y],'kfold',KF,'options',opt);
        n
        tr
    end
end
vals = mean(vals,4);
%% display hold out over factors
close all
ME = squeeze(mean(vals(:,1,:),1));
sME = squeeze(std(vals(:,1,:),1));
ME = ME.^.5;
sME = sME.^.5;
%sME = sME*(ITR*KF)^-.5;
figure;
plot(ME);
hold on
plot(ME+sME,'r--')
hold on
plot(ME-sME,'r--')
RHO = squeeze(mean(vals(:,2,:),1)); 
sRHO = squeeze(std(vals(:,2,:),1));
sRHO = sRHO*(ITR*KF)^-.5;
figure
plot(RHO);
hold on
plot(RHO+sRHO,'r--')
hold on
plot(RHO-sRHO,'r--')

[XL,YL,XS,YS,BETA,PCTVAR,MSE,st] = plsregress(X,Y,10);
Ya = [ones(size(X,1),1) X]*BETA;
Ya = PCA_BKPROJ(Ya,gE,gU);
[R P] = corr(gS,Ya);
figure;plot(diag(R));
title(['Correlation--' num2str(R) '--Significance--' num2str(P)]);
ylabel('Yield Approximate','fontsize',10);
xlabel('Yield Actual','fontsize',10);
figure;
plot(BETA(2:end));
%% column wise
Ya = [];
for e = 1:size(Y,2)
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,st] = plsregress(X,Y(:,e),10);
    Ya(:,e) = [ones(size(X,1),1) X]*BETA;
end
Ya = PCA_BKPROJ(Ya,gE,gU);
figure;
[R P] = corr(gS,Ya);
figure;plot(diag(R));
hold on;