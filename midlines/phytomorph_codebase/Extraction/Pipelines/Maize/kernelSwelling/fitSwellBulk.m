function [p ep er] = fitSwellBulk(d,ft,TAU)
    if nargin == 2
        TAU = 1:size(d,1);
    end
    numberOfVariables = 2;
    options = optimset('MaxFunEvals',3000*numberOfVariables,'TolFun',10^-6);
    for tr = 1:size(d,2)        
        %[p(tr,:) er(tr)] = fminsearch(@(X)mySwellFit(d(:,tr),X),[10^4 .01]);
        FINIT = mean(d(:,tr)');
        [p(tr,:) er(tr) ef(tr) output{tr}] = fminsearch(@(X)mySwellFit(d(:,tr)',X),[FINIT .01],options);
        ep(tr,:) = func(p(tr,1),p(tr,2),0:(ft-1));
    end
end

%{
    % test with Jeff's data
    swellD = readtext(['/home/nate/Downloads/Kernel_LogFits_byGenotype.csv']);
    swellD = readtext(['/home/nate/Downloads/Kernels_PercentIncrease_byGenotype.csv']);
    G = swellD(2:end,1);
    swellD = cell2mat(swellD(2:end,2:end));
    swellD = swellD';
    TAU = 20:size(swellD,1);
    [p ep er] = fitSwellBulk(swellD,size(swellD,1),TAU);
    close all
    figure;
    plot(swellD,'k')
    hold on
    plot(ep','r')
    

    UQ = unique(G);
    figure;
    hold on

    for u = 1:numel(UQ)
        uidx = strcmp(G,UQ{u});
        sD = ep(uidx,:);
        U = mean(sD,1);
        SE = std(sD,1,1).*size(sD,1)^-.5;
        errorbar(U,SE);
    end
    legend(UQ)
    
    
    TTPV = table;
    TTh = table;
    tp = p;
    tp = bp;
    for u1 = 1:numel(UQ)
        uidx1 = strcmp(G,UQ{u1});
        for u2 = 1:numel(UQ)
            uidx2 = strcmp(G,UQ{u2});
            
            d1 = tp(uidx1,1);
            TF1 = isoutlier(d1);

            d2 = tp(uidx2,1);
            TF2 = isoutlier(d2);

            %d1(TF1) = [];
            %d2(TF2) = [];

            [h,pv] = ttest2(d1,d2);
            TTPV{UQ{u1},UQ{u2}} = pv;
            TTh{UQ{u1},UQ{u2}} = h;
        end
    end


    TTPV2 = table;
    TTh2 = table;
    for u1 = 1:numel(UQ)
        uidx1 = strcmp(G,UQ{u1});
        for u2 = 1:numel(UQ)
            uidx2 = strcmp(G,UQ{u2});
            
            d1 = tp(uidx1,2);
            TF1 = isoutlier(d1);

            d2 = tp(uidx2,2);
            TF2 = isoutlier(d2);

            %d1(TF1) = [];
            %d2(TF2) = [];
            TTPV2{UQ{u1},UQ{u2}} = pv;
            TTh2{UQ{u1},UQ{u2}} = h;
        end
    end

    %%
    rmfp = [];
    fV = {};
    parfor u = 1:numel(UQ)
        tic
        uidx = strcmp(G,UQ{u});
        subD = swellD(:,uidx)';
        [fixed{u},random{u},stats{u}] = fitSwellCurvePerGroup(subD);
        mp = bsxfun(@plus,random{u},fixed{u})';
        for i = 1:size(mp,1)
            fV{u}(:,i) = func(mp(i,1),mp(i,2),1:size(subD,2));
        end
        toc
    end

    cnt = 1;
    bp = [];
    for u = 1:numel(UQ)
        mp = bsxfun(@plus,random{u},fixed{u})';
        bp = [bp;mp];
        for i = 1:size(fV{u},2)
            betterF(cnt,:) = fV{u}(:,i)';
            cnt = cnt + 1;
        end
    end

    outl1 = isoutlier(bp(:,1));
    outl2 = isoutlier(bp(:,2));
    rm = find(outl1 | outl2);


    outL = swellD(:,rm);

    G(rm) = [];
    swellD(:,rm) = [];
    figure;
    plot(mean(swellD,2));
    hold on
    plot(mean(outL,2))

    


    rmfp = [];
    fV = {};
    parfor u = 1:numel(UQ)
        tic
        uidx = strcmp(G,UQ{u});
        subD = swellD(:,uidx)';
        [fixed{u},random{u},stats{u}] = fitSwellCurvePerGroup(subD);
        mp = bsxfun(@plus,random{u},fixed{u})';
        for i = 1:size(mp,1)
            fV{u}(:,i) = func(mp(i,1),mp(i,2),1:size(subD,2));
        end
        toc
    end

    cnt = 1;
    bp = [];
    for u = 1:numel(UQ)
        mp = bsxfun(@plus,random{u},fixed{u})';
        bp = [bp;mp];
        for i = 1:size(fV{u},2)
            betterF(cnt,:) = fV{u}(:,i)';
            cnt = cnt + 1;
        end
    end

    
   




%}