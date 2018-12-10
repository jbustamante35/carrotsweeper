D = readtext('/home/nate/Downloads/cotdataRaw.csv');
D = readtext('/home/nate/Downloads/cot_train_flg22.csv');
HEAD = D(1:2,:);
X = cell2mat(D(4:end,:));
Y = cell2mat(D(3,:));
[cotS cotC cotU cotE cotL cotERR cotLAM] = PCA_FIT_FULL_T(X,3);
%%
pD = readtext('/home/nate/Downloads/cotdataRaw_predict.csv');
D = readtext('/home/nate/Downloads/cot_predict_flg22.csv');
pHEAD = pD(1:2,:);
pX = cell2mat(pD(4:end,:));
pY = cell2mat(pD(3,:));
%% combine Header
for e = 1:size(pHEAD,2)
    cHEAD{e} = [pHEAD{1,e} pHEAD{2,e}];
end

for u = 1:numel(cHEAD)
    cHEAD{u} = strrep(strrep(cHEAD{u},' ',''),'/','');
end
%%
close all
UQ = unique(Y);
CL = {'r.' 'b.'}
for u = 1:numel(UQ)
    fidx = find(Y==UQ(u));
    plot3(cotC(1,fidx),cotC(2,fidx),cotC(3,fidx),CL{u});
    hold on
end
%%
close all
lambda = myLDA(cotC',Y);
Flambda = PCA_BKPROJ_T(lambda,cotE,cotU);
scores = lambda'*cotC;
for u = 1:numel(UQ)
    ksdensity(scores(Y==UQ(u)));
    hold on
end
%%
close all
plot(Flambda - cotU)
%%
[B,FitInfo] = lasso(X',Y,'CV',5);
%%
aX = pX;
aY = pY;
aHEAD = cHEAD;
close all
lassoPlot(B,FitInfo,'PlotType','CV');
scores = aX'*B(:,FitInfo.IndexMinMSE);
CL = {'k' ,'r'};
LABEL{1} = {'non-leasion'};
LABEL{2} = {'leasion'};
subGRPS = 2;

UQ = unique(Y);
for u = 1:numel(UQ)
    fidx = find(aY==UQ(u));
    sig = scores(aY==UQ(u));
    [f,xi] = ksdensity(sig);
    plot(xi,f,CL{u})
    kidx = kmeans(sig,subGRPS);
    
    
    for k = 1:subGRPS
        fprintf('***********************')
        fprintf(['group:' num2str(u) '\n'])
        hGRP{u,k}.grpMean = mean(sig(kidx==k));
        hGRP{u,k}.index = fidx(kidx==k);
        [hGRP{u,k}.uniqueHeader ia ic] = unique(aHEAD(:,fidx(kidx==k)));
        mu = unique(ic);
        for uu = 1:numel(mu)
            hGRP{u,k}.NUM(uu) = sum(ic == mu(uu));
        end
        hGRP{u,k}.TYPE = LABEL{u};
    end
    
    %{
    fprintf('***********************')
    mean(sig(kidx==2))
    fidx(kidx==2)
    aHEAD(:,fidx(kidx==2))
    fprintf('***********************')
    %}
    hold on
end
%%
ML = unique(aHEAD);
SHEET = table;
cnt = 1;
for u = 1:numel(ML)
    ML{u} = strrep(strrep(ML{u},' ',''),'/','');
end

for u = 1:numel(ML)
    cnt = 1;
    for h1 = 1:size(hGRP,1)
        for h2 = 1:size(hGRP,2)
            SHEET{cnt,ML{u}} = 0;
            cnt = cnt + 1;
        end
    end
end
%%
ROWN = {};
for u = 1:numel(ML)
    cnt = 1;
    for h1 = 1:size(hGRP,1)
        for h2 = 1:size(hGRP,2)
            
            for h = 1:numel(hGRP{h1,h2}.uniqueHeader)
                if strcmp(ML{u},hGRP{h1,h2}.uniqueHeader{h})
                    SHEET{cnt,strrep(ML{u},' ','')} = SHEET{cnt,strrep(ML{u},' ','')} + hGRP{h1,h2}.NUM(h);
                end
            end
             ROWN{cnt} = [hGRP{h1,h2}.TYPE{1} '--' num2str(hGRP{h1,h2}.grpMean)];
            cnt = cnt + 1;
            
           
        end
    end
end
SHEET.Row = ROWN;

%%
close all
plot(B(:,FitInfo.IndexMinMSE))