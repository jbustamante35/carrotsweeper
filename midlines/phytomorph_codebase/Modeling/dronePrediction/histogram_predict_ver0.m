% NOTES
% check to see if plots swtich postions in the sheets
%%
FilePath = '/mnt/spaldingdata/Drone Imagery/Arlington_2017/NDVI_For_Nathan/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
pY = readtext(FileList{end});
sT = cell2mat(pY(2:end,2));
tT = cell2mat(pY(2:end,3));
%% 
GY = [];
for s = 1:(numel(FileList)-1)
    d = readtext(FileList{s});
    %% correct the plot name of border
    flag = 0;
    for e = 2:(size(d,2))
        if ischar(d{1,e})
            if flag == 0
                d{1,e} = 125;
                flag = 1;
            else
                d{1,e} = 225;
            end
        end

    end
    %%
    H = cell2mat(d(1,2:end));
    H
    nH = 1:numel(H);
    X = cell2mat(d(2:end,1));
    Y = cell2mat(d(2:end,2:end));
    GY = cat(3,GY,Y);
    fprintf(['done reading sheet:' num2str(e) ':' str2num(numel(FileList)-1) '\n'])
end
%% normalize
for p = 1:size(GY,2)
    for t = 1:size(GY,3)
        nGY(:,p,t) = GY(:,p,t)*sum(GY(:,p,t))^-1;
    end
end
%% smoooth and normalize
N = 101;
for p = 1:size(GY,2)
    for t = 1:size(GY,3)
        sig = GY(:,p,t);
        sig = imfilter(sig,ones(N,1)/N,'replicate');
        snGY(:,p,t) = sig*sum(sig)^-1;
    end
end
%% look at single plot
plotNum = 10;
toPlot = snGY;
close all
figure
plot(squeeze(toPlot(:,plotNum,:)));
figure
for t = 1:size(toPlot,3)
    plot3(X,toPlot(:,plotNum,t),t*ones(size(toPlot,1),1))
    hold on
end
%% decompose
sz = size(snGY);
tmpD = reshape(snGY,[prod(sz(1)) prod(sz(2:3))]);
[S C U E L ERR LAM] = PCA_FIT_FULL_T(tmpD,size(tmpD,1));
[S C U E L ERR LAM] = PCA_FIT_FULL_T(tmpD,4);
dD = reshape(C,[size(C,1) sz(2:3)]);
rX = permute(dD,[2 1 3]);
szR = size(rX);
rX = reshape(rX,[szR(1) prod(szR(2:3))]);
%% sweep
tmpC = mean(C,2);
for c = 1:size(tmpC,1)
    L = linspace(min(C(c,:)),max(C(c,:)),5);
    for l = 1:numel(L)
        tC = tmpC;
        tC(c) = L(l);
        M = PCA_BKPROJ_T(tC,E,U);
        plot(X,M);
        hold on
    end
    waitforbuttonpress
    hold off
end
%% regression
close all
h1 = figure;
h2 = figure;
rX = rX;
%rX = profileTOPX1000;
rY = tT;
COR = [];
Yp = [];
for l = 1:10
    for e = 1:size(rX,1)
        idx = setdiff(1:size(rX,1),e);
        
        subX = rX(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        %size(R);
        subY = rY(idx);
        
        %{
        %subX = [Co(idx,:) LV(idx,1)];
        subX = [Co(idx,:)];
        net = feedforwardnet([l 5]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        Yp(l,e) = sim(net,Co(e,:)');
        %}
        
        
        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
        %Yp = [ones(size(X,1),1) X(e,:)]*beta;
        %Yp(e) = [1 X(e,:)]*beta;
        %Yp(l,e) = [1 X(e,:) LV(e,1)]*beta;
        Yp(l,e) = [1 rX(e,:)]*beta;
        %{
        figure(h2)
        plot(rY(e),Yp(l,e),'.')
        drawnow
        hold on
        %}
    end
    COR(l) = corr(rY(:),Yp(l,:)','type','Pearson');
    %l
    figure(h1)
    plot(COR)
    drawnow
end
%% tassel time ridge
figure;
centerX = bsxfun(@minus,rX,mean(rX,1));
uY = mean(tT);
centerY = tT - uY;
close all
options = statset('UseParallel',true);
[B FitInfo] = lassoglm(centerX,centerY,'normal','CV',10,'Alpha',.0001,'Options',options);
lassoPlot(B,FitInfo,'plottype','CV');
find(B(:,FitInfo.Index1SE))
preY = centerX*B(:,FitInfo.Index1SE) + uY;
figure
plot(tT,preY,'.')
title(['Tassel Time-Ridge' num2str(corr(preY,tT))]);

%% tassel time lasso
close all
figure;
centerX = bsxfun(@minus,rX,mean(rX,1));
uY = mean(tT);
centerY = tT - uY;
options = statset('UseParallel',true);
[B FitInfo] = lassoglm(centerX,centerY,'normal','CV',10,'Alpha',1,'Options',options);
lassoPlot(B,FitInfo,'plottype','CV');
find(B(:,FitInfo.Index1SE))
preY = centerX*B(:,FitInfo.Index1SE) + uY;
figure
plot(tT,preY,'.')
title(['Tassel Time-Lasso' num2str(corr(preY,tT))]);


%% tassel time elastic net
close all
figure;
centerX = bsxfun(@minus,rX,mean(rX,1));
uY = mean(tT);
centerY = tT - uY;
options = statset('UseParallel',true);
[B FitInfo] = lassoglm(centerX,centerY,'normal','CV',10,'Alpha',.5,'Options',options);
lassoPlot(B,FitInfo,'plottype','CV');
find(B(:,FitInfo.Index1SE))
preY = centerX*B(:,FitInfo.Index1SE) + uY;
figure
plot(tT,preY,'.')
title(['Tassel Time-elastic net' num2str(corr(preY,tT))]);

%% silk time ridge
centerX = bsxfun(@minus,rX,mean(rX,1));
uY = mean(sT);
centerY = sT - uY;
close all
options = statset('UseParallel',true);
[B FitInfo] = lassoglm(centerX,centerY,'normal','CV',10,'Alpha',.0001,'Options',options);
lassoPlot(B,FitInfo,'plottype','CV');
find(B(:,FitInfo.Index1SE))
preY = centerX*B(:,FitInfo.Index1SE) + uY;
figure
plot(sT,preY,'.')
title(['Tassel Time-Ridge' num2str(corr(preY,sT))]);
%{
%% silk time lasso
centerX = bsxfun(@minus,rX,mean(rX,1));
uY = mean(sT);
centerY = sT - uY;
close all
options = statset('UseParallel',true);
[B FitInfo] = lassoglm(centerX,centerY,'normal','CV',10,'Alpha',1,'Options',options);
lassoPlot(B,FitInfo,'plottype','CV');
find(B(:,FitInfo.Index1SE))
preY = centerX*B(:,FitInfo.Index1SE) + uY;
figure
plot(sT,preY,'.')
title(['Tassel Time-Lasso' num2str(corr(preY,sT))]);
%}
%% tassel time elastic net
close all
figure;

%tmpX = rX(:,(end-100):end,1);
tmpX = rX(:,:);
%tmpX = sort(tmpX,2);




kidx = find(mean(tmpX(:,1:10),2) < 30);
tmpX = tmpX(kidx,:);
tmpY = rY(kidx);
centerX = bsxfun(@minus,tmpX,mean(tmpX,1));
uY = mean(tmpY);
centerY = tmpY - uY;
options = statset('UseParallel',true);
[B FitInfo] = lassoglm(centerX,centerY,'normal','CV',10,'Alpha',1,'Options',options);
lassoPlot(B,FitInfo,'plottype','CV');
find(B(:,FitInfo.Index1SE))
preY = centerX*B(:,FitInfo.Index1SE) + uY;
figure
plot(tmpY,preY,'.')
title(['Tassel Time-elastic net' num2str(corr(preY,tmpY))]);


%%
%%
close all
plot(mean(tmpX(kidx,1:10),2),rY(kidx),'.')
corr(mean(tmpX(kidx,1:10),2),rY(kidx))
%%

rm = find(mean(tmpX(:,1:10),2) > 30);



tmpX(rm,:) = [];
rY(rm) = [];



















