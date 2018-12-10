D{1} = csvread('/home/nate/Downloads/esb1.csv');
D{2} = csvread('/home/nate/Downloads/WTCol.csv');
D{1} = D{1}(:,2:end);
D{2} = D{2}(:,2:end);
[h] = ttest2(D{1},D{2},[],[],[],2);
W = {};
U = {};
S = {};
G = [];
DM = [];
figure;
hold all
scales = [1:20];
for e = 1:numel(D)
    U{e} = mean(D{e},2);
    S{e} = std(D{e},1,2).*size(D{e},2).^-.5;
    errorbar(U{e},S{e})
    LEG{e} = num2str(e);
    for k = 1:size(D{e},2)
        W{e}(:,:,k) = cwt(D{e}(:,k),scales,'gaus1')';
    end
    DM = [DM D{e}];
    G = [G e*ones(1,size(D{e},2))];
end
plot(h,'r')
legend(LEG);
for w = 1:numel(scales)
    H(:,w) = ttest2(squeeze(W{1}(:,w,:)),squeeze(W{2}(:,w,:)),[],[],[],2);
end
plot(any(H,2),'c');
%%
csvwrite('/mnt/spaldingdata/nate/ttest.csv',h)
%% workup with splitting vector
close all
[lambda] = myLDA(DM',G');
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(DM',2);
[plambda] = myLDA(sC,G');
plambdaF = PCA_BKPROJ(plambda',sE,sU);
plot(plambdaF-sU);
sc = plambdaF*bsxfun(@minus,DM,sU');
sc = (plambdaF-sU)*bsxfun(@minus,sS',sU');
%sc = plambdaF*sS'
[hl pl] = ttest2(sc(G==1),sc(G==2));
sc2 = plambda'*sC';
[hl2 pl2] = ttest2(sc2(G==1),sc2(G==2));
figure;
UQ = unique(G);
CL = {'r.','b.'};
for u = 1:numel(UQ)
    plot(sC(G==UQ(u),1),sC(G==UQ(u),2),CL{u});
    hold on
end
quiver(0,0,plambda(1),plambda(2),30,'g')
%%
close all
[mlambda] = mynLDA(sC,G',1,2);
mplambdaF = PCA_BKPROJ(mlambda',sE,sU);
plot(bsxfun(@minus,mplambdaF,sU)');