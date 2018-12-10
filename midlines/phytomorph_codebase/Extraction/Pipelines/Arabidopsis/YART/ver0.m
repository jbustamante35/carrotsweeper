FilePath = '/mnt/spaldingdata/nate/mirror_images/forRichard/';
FileList = {};
FileExt = {'tiff','TIF','tif','JPG'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% view movie
for e = 1:numel(SET)
    for f = 1:numel(SET{e})
        I = imread(SET{e}{f});
        imshow(I,[]);
        drawnow
    end
end
%% read stack
for e = 1:numel(SET)
    for f = 1:numel(SET{e})
        ST(:,:,:,f) = imread(SET{e}{f});
        f
    end
end
%%
myTraceLow(ST(:,:,1),1)
%%
para.scales.value=5;
para.resize.value=1;
K = surKur(ST(:,:,1),para);
%% threshold the stack on curvature
tmp = bindVec(K(:,:,1));
B = tmp > graythresh(tmp);
[x1 x2 V] = impixel(ST(:,:,1));
SKEL = bwmorph(B,'skel',inf);

[x y] = find(SKEL);
DP = [x y]';
T = Radjacency(DP,10);
for e = 1:size(x1,1)
    delta = (x - x2(e)).^2 + (y - x1(e)).^2;
    [~,idx(e)] = min(delta);
end
[path , pathcost]  = dijkstra(T , idx(1) , idx(2));
gamma = DP(:,path);
%%

imshow(ST(:,:,1),[])
hold on
plot(gamma(2,:),gamma(1,:),'r')

%%
imshow(cat(3,ST(:,:,1,1),K));
%% find corners
[cim, r1, c1] = harris(ST(:,:,1,1),5, .05, 20,1);
%% align
[optimizer, metric] = imregconfig('monomodal')
for f = 1%:size(ST,3)-1
    fixed = rgb2gray(ST(:,:,:,f));
    move = rgb2gray(ST(:,:,:,f+1));
    [F] = stillMain(fixed,move,1);
    
    %reg = imregister(move, fixed, 'affine', optimizer, metric,'DisplayOptimization',1); 
end






