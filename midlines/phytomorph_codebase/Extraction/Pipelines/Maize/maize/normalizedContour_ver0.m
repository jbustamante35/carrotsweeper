C = [];
tmpK = [];
for e = 1:numel(moleculeStore.store)
        try
        [~,~,tmpK(:,:,e),~] = moleculeStore.store{e}.toTensorForm();
        T = tmpK(:,end,e) - tmpK(:,1,e);
        U = mean([tmpK(:,end,e) tmpK(:,1,e)],2);
        L = norm(T);
        T =  T / L;
        N = [T(2);-T(1)];
        T = T / L;
        T = T * 150;
        E = [T N];
        tmpC = PCA_REPROJ(tmpK(:,:,e)',E,U');
        hold on;
        plot(tmpC(:,1),tmpC(:,2))
        tmpC = tmpC';
        C = cat(3,C,tmpC);
        tmpC(1,:) = -tmpC(1,:);
        tmpC = fliplr(tmpC);
        C = cat(3,C,tmpC);
        catch
        end
end
%% spline
figure;
clear kC;
for e = 1:size(C,3)
    fn{e} = spap2(10,3,linspace(0,1,size(C,2))',C(:,:,e));
    sf = fnval(fn{e},linspace(0,1,size(C,2))');
    kC(e,:) = fn{e}.coefs(:);
    plot(sf(1,:),sf(2,:))
    hold on;
    plot(C(1,:,e),C(2,:,e),'r');
    hold off
    drawnow
    pause(.1);
end
%%
figure;
hold on;
for e = 1:size(C,3)
    plot(C(1,:,e),C(2,:,e))
end
%%
sz = size(C);
C = reshape(C,[size(C,1)*size(C,2) size(C,3)]);
[kS kC kU kE kL kERR kLAM] = PCA_FIT_FULL(C',6);
kS = reshape(kS',sz);
C = reshape(C,sz);
%%
for e = 1:size(kS,3)
    plot(kS(1,:,e),kS(2,:,e),'r')
end
%%
kU = reshape(kU,sz(1:2));
figure;
plot(kU(1,:),kU(2,:));
kE = reshape(kE,[sz(1:2) size(kE,2)]);
%hold on
for comp = 1:size(kE,3)
   
    quiver(kU(1,:),kU(2,:),kE(1,:,comp),kE(2,:,comp),1)
    waitforbuttonpress
end
%%
%% get error contour
[~,~,k,~] = curves{e}(c).toTensorForm();
T = k(:,end) - k(:,1);
U = mean([k(:,end) k(:,1)],2);
L = norm(T);
T =  T / L;
N = [T(2);-T(1)];
T = T / L;
E = [T N];
k = PCA_REPROJ(k',E,U');
k = k';

szk = size(k);
rk = k;
k = reshape(k,[1 numel(k)]);
kC = PCA_REPROJ(k,kE,kU);
k = PCA_BKPROJ(kC,kE,kU);
k = reshape(k',szk);
figure;
plot(k(1,:),k(2,:));
hold on;
plot(rk(1,:),rk(2,:),'r')
%%
% generate the bounds of the data
close all
mag = 1;
uMIN = mag*min(kC,[],1);
uMAX = mag*max(kC,[],1);            
b = [[uMIN;uMAX]];
POP = 50;
GEN = 1000;
options = psooptimset('PopInitRange',b,'PopulationSize',POP,'Display','iter','Generations',GEN,'Vectorized','on','CognitiveAttraction',2,'SocialAttraction',.1);    
% run the optimization
percent = .85;
[~,~,k,~] = curves{end}(c).toTensorForm();
T = k(:,end) - k(:,1);
U = mean([k(:,end) k(:,1)],2);
L = norm(T);
T =  T / L;
N = [T(2);-T(1)];
T = T / L;
T = T * 150;P = intper
E = [T N];
k = PCA_REPROJ(k',E,U');
targetCurve = k';
func = @(X)PCDs(X,kU,kE,targetCurve,percent,10,fn{1});
[mo,fval,exitflag,output,population,scores] = pso(func,size(b,2),[],[],[],[],[],[],[],options);
%sourceCurve = PCA_BKPROJ(mo,kE,kU);
%sourceCurve = reshape(sourceCurve,size(targetCurve));
spline = fn{1};
spline.coefs = reshape(mo,size(spline.coefs));
sourceCurve = fnval(spline,linspace(0,1,size(targetCurve,2)));
close all
plot(targetCurve(1,:),targetCurve(2,:),'r');
hold on
plot(sourceCurve(1,:),sourceCurve(2,:))
%%
close all
[mo,fval,exitflag] = fminsearch(func,zeros(1,size(kC,2)));
%%
close all
[mo,fval,exitflag] = fminunc(func,zeros(1,size(kC,2)));
hold on
plot(sourceCurve(1,:),sourceCurve(2,:))