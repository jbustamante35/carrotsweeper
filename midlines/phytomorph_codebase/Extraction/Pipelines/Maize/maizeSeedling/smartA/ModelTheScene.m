FilePath = '/mnt/tetra/nate/sceneModel/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
toLoad = 140;
toLoad = numel(FileList);

k=1;
e = 1;
a = load(FileList{e},'RESIZED_sceneStruct');
tmp = a.RESIZED_sceneStruct.conetainerBox(k).IMAGE;
tmp = imresize(tmp,.15);
newSZ = size(tmp);

for e = 1:toLoad
    MS{e} = zeros([newSZ 3]);
    e
end

parfor e = 1:toLoad
    try
        tic
        a = load(FileList{e},'RESIZED_sceneStruct');
        n = numel(a.RESIZED_sceneStruct.conetainerBox);
        %MS{e} = zeros([newSZ 3]);
        
        for k = 1:n
            tmp = a.RESIZED_sceneStruct.conetainerBox(k).IMAGE;
            CEN{e}(k,:) = a.RESIZED_sceneStruct.conetainerBox(k).Centroid;
            LOC{e}(k,:) = a.RESIZED_sceneStruct.conetainerBox(k).coneTainerLocation;
            %containerStack(:,:,:,cnt) = imresize(tmp,.5);
            MS{e}(:,:,:,k) = imresize(tmp,.15);
            %imshow(containerStack(:,:,:,cnt),[]);
            %drawnow
            %cnt = cnt + 1;
        end
        toc
        e
    catch ME
        ME
    end
end
%%
mssz = size(MS{1});
containerStack = zeros([mssz(1:3),toLoad*3]);
containerIndex = zeros(toLoad*3,1);
containerLocation = zeros(toLoad*3,2);
cnt = 1;
for e = 1:numel(MS)
    for k = 1:size(MS{e},4)
        if size(MS{e},4) == size(LOC{e},1) & size(MS{e},4) == size(CEN{e},1)
            containerStack(:,:,:,cnt) = MS{e}(:,:,:,k);
            containerIndex(cnt,:) = LOC{e}(k);
            containerLocation(cnt,:) = CEN{e}(k,:);
            cnt = cnt + 1;
        end
    end
end
if cnt <= size(containerStack,4)
    containerStack(:,:,:,cnt:end) = [];
    containerIndex(cnt:end) = [];
    containerLocation(cnt:end,:) = [];
end
%% mean has label
close all
sz = size(containerStack);
d = reshape(containerStack,[prod(sz(1:3)) sz(4)]);
%% PCA data
kp = 1:size(d,2);
kp = find(cERR < 13);
%%
[cU,cE,cL] = PCA_FIT_FULL_Tws(d(:,kp(1:200)),5);
cC = PCA_REPROJ_T(d(:,kp),cE,cU);
cSIM = PCA_BKPROJ_T(cC,cE,cU);
cERR = sum((d(:,kp) - cSIM).^2,1).^.5;
%kp = find(cERR < 20);

cSIM = reshape(cSIM,[sz(1:3) numel(kp)]);
uSIM = reshape(cU,sz(1:3));
imshow(uSIM,[]);
%%
close all
figure
for e = 1:100
    imshow([cSIM(:,:,:,e) containerStack(:,:,:,kp(e))],[])
    drawnow
    pause(.2)
end
%% cca
X = d';
X = cC';
Y = containerLocation;
[A,B,r,U,V,stats] = canoncorr(X,Y);

%% sweep location cca
close all

[s1 s2] = ndgrid(linspace(min(containerLocation(:,1)),max(containerLocation(:,1)),100),linspace(mean(containerLocation(:,2)),mean(containerLocation(:,2)),1));

for e = 1:size(s1,1)
    co = [s1(e) s2(e) 1];
    coS = co*b;
    coS = PCA_BKPROJ_T(coS',cE,cU);
    coS = reshape(coS,sz(1:3));
    imshow(coS,[]);
    drawnow
end

%% regress location from data
close all
X = d';
Y = containerLocation;
%B = lasso(X,Y);
[Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(X,Y,15);
preLoc = [ones(size(d,2),1) d']*beta;
figure;
plot(preLoc(:,1),containerLocation(:,1),'.');
figure;
plot(preLoc(:,2),containerLocation(:,2),'.');

%% regress container based on location
close all
b = [];
preY = [];
for pc = 1:size(cC,1)
    Y = cC(pc,:)';
    X = [containerLocation ones(size(containerLocation,1),1)];
    [b(:,pc),bint,r,rint,stats] = regress(Y,X);
    
    preY(:,pc) = X*b(:,pc);
    
    
    
    plot(X(:,1),Y,'b.')
    hold on
    plot(X(:,1),preY(:,pc),'c.')
    hold on
    plot(X(:,2),Y,'r.')
    plot(X(:,2),preY(:,pc),'m.')
    hold off
    pause(.2)
end

lSIM = PCA_BKPROJ_T(preY',cE,cU);
lSIM = reshape(lSIM,sz);
%% sweep location
close all
[s1 s2] = ndgrid(linspace(min(containerLocation(:,1)),max(containerLocation(:,1)),100),linspace(mean(containerLocation(:,2)),mean(containerLocation(:,2)),1));
for e = 1:size(s1,1)
    co = [s1(e) s2(e) 1];
    coS = co*b;
    coS = PCA_BKPROJ_T(coS',cE,cU);
    coS = reshape(coS,sz(1:3));
    imshow(coS,[]);
    drawnow
end








