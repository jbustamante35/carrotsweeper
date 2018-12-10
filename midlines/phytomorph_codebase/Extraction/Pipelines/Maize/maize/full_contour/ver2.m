FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%%
set = randi(numel(SET),1);
I = double(imread(SET{set}{1}));
[curveBank] = getLevelContours(I,30);
[curveBank] = selectClosedCurves(curveBank);
rmidx = [curveBank.length] < 200 | [curveBank.length] > 1000;
curveBank(rmidx) = [];
nSZ = [20 20];
[n1 n2] = ndgrid(linspace(-nSZ(1),nSZ(1),2*nSZ(1)),linspace(-nSZ(2),nSZ(2),2*nSZ(2)));
sz = size(n1);
NH = [n1(:) n2(:)]';
tic
[curveBank CB] = sampleCurveBank(I,curveBank,NH,sz);
toc
%% populate curve segments
for e = 1:numel(CB)
    for e1 = 1:numel(CN{e})
        D(cnt,:) = CB{e}{e1}.distance
    end
end
%%
close all
imshow(I,[])
hold on
for e = 1:numel(curveBank)
    plot(curveBank(e).data(1,:),curveBank(e).data(2,:),'r');
    drawnow
end
%% 
SZ = 51;
SNIP = 10;
BANK = zeros(2*SZ,100000);
IDX = zeros(2*SZ,100000);
cnt = 1;
%h1 =  figure;
%h2 =  figure;

for e = 1:numel(curveBank)
    tmp = [curveBank(e).data];
    tmpSZ = size(tmp,2);
    %{
    figure(h1);
    plot(tmp(1,:),tmp(2,:),'r')
    hold on
    %}
    for e1 = 1:size(tmp,2)
        sig = tmp(:,1:SZ)';
        fn = spap2(1,3,[1:size(sig,1)]',sig');
        fn2 = fnder(fn,1);
        
        tvec = fnval(fn2,(SZ-1)/2);
        tvec = tvec/norm(tvec);
        nvec = [tvec(2);-tvec(1)];
        
        E = [tvec nvec];
        U = fnval(fn,(SZ-1)/2);
        C = PCA_REPROJ(sig,E,U');
        BANK(:,cnt) = C(:);
        IDX(:,cnt) = e*ones(size(C(:)));
        cnt = cnt + 1;
        
        
        %BANK = [BANK C(:)];
        %{
        figure(h1);
        plot(sig(:,1),sig(:,2),'r');
        hold on
        quiver(sig((SZ-1)/2,1),sig((SZ-1)/2,2),nvec(1),nvec(2),10,'b')
        quiver(sig((SZ-1)/2,1),sig((SZ-1)/2,2),tvec(1),tvec(2),10,'g')        
        drawnow
        
        
        
        figure(h2);
        plot(C(:,1),C(:,2),'r')
        hold on
        drawnow
        %}
        
        tmp = circshift(tmp,[0 1]);
        
    end
    e
end
%% remove empty from bank
ridx = find(all(BANK == 0,1));
BANK(:,ridx) = [];
IDX(:,ridx) = [];
IDX = mean(IDX,1);
%%
[S C U E L ERR LAM] = PCA_FIT_FULL(BANK',3);
UQ = unique(IDX)
h1 = figure;
h2 = figure;
for u = 1:numel(UQ)
    fidx = find(UQ(u)==IDX);
    sig = C(fidx,:);
    
    figure(h2);
    plot3(sig(:,1),sig(:,2),sig(:,3),'r')
    
    
    figure(h1)
    plot(curveBank(UQ(u)).data(1,:),curveBank(UQ(u)).data(2,:));
    
    pause(.5);
end

plot3(sig(:,1),sig(:,2),sig(:,3),'b');
%%
D = pdist(C);
D = squareform(D);
