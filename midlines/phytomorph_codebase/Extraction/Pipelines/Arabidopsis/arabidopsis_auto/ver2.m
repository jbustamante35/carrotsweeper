%%
FilePath = '/mnt/spaldingimages/nate/whole_Candace_Gravitropism_arabidopsis/';
FilePath = '/mnt/scratch1/myLinks/arabidopsisGravitropism/baseLine/';
FilePath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/toProcess/';
FilePath = '/mnt/scratch1/myLinks/arabidopsisGravitropism/candace_qtl/';
FilePath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/toProcess/straight/RafaelStraight14Aug13/lip5-1_Cam2/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% sort sets
for e = 1:numel(SET)
    NAME = [];
    for i = 1:numel(SET{e})
        [pth,nm,ext] = fileparts(SET{e}{i});
        NAME(i) = str2num(nm);
    end
    [J sidx] = sort(NAME);
    SET{e} = SET{e}(sidx);
    e
end
%% try again
oPath = '/mnt/scratch5/arabidopsis_model_13.09.09/';
oPath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/new_return/';
mkdir(oPath);
for e = 1:numel(SET)
    tm = clock;
    isolateRoots_overStack(SET{e},oPath,1,3,20,20,20,1);
    etm = etime(clock,tm);
    etm
end
%%

%%
for e = 1:numel(out)
    angle(:,e) = measureStruct(out{e},20,[301 1]);
end
%% watch stack
for e = 1:numel(SET)
    out{e} = isolateRoots_overStack(SET{e},'/mnt/scratch5/arabidopsis_model_13.09.06/',0);
end
%% view stack
for e = 1:10%numel(SET)
    viewStruct(out{e},SET{e});
end
%% first images
for e = 1:numel(SET)
    isolateRoots(SET{e}{1},1);
end
%% view
for s = 1:numel(SET)
    for e = 1:numel(SET{s})
        I = imread(SET{s}{e});
        imshow(I,[])
        drawnow
    end
end
%% diff
for s = 1:numel(SET)    
    I1 = imread(SET{s}{1});
    Ie = imread(SET{s}{end});
    D = Ie - I1;
    imshow(D,[])
    drawnow
end
%%
fn = surfaceKurvatureFeatureMap();
t = imread(SET{1}{end});
h = fspecial('gaussian', 11,6);
t = imfilter(t,h);
sz = size(t);
to = t;
t = imresize(t,.25);
para.scales.value = 1;
para.resize.value = 1;
[K,Vo] = surKur(t,para);

V1 = imresize(Vo(:,:,1),sz);
V2 = imresize(Vo(:,:,2),sz);

W1 = imresize(Vo(:,:,3),sz);
W2 = imresize(Vo(:,:,4),sz);
%%
imshow(to,[])
hold on
quiver(V1,V2,'b')
quiver(W1,W2,'r')
%%
%[c r V] = impixel(to);
close all
c = (1:1:size(to,2))';
%c = linspace(200,size(to,2)-200,2000)';
r = (size(to,1)-100*ones(to,size(c,1)))';
%P = [c,r];
%imshow(to);hold on
for p = 2:size(c,1)-1
   P(p).path = [c(p-1:p+1),r(p-1:p+1)'];
   [P(p).path P(p).dist] = traceKpath(P(p).path,W1,W2);
   p
end
%%
close all
imshow(to);hold on
for p = 1:numel(P)
    try
        if P(p).dist(end) < 7
            plot(squeeze(P(p).path(2,1,:)),squeeze(P(p).path(2,2,:)));
        end
    catch ME
    end
end
     %{   
    [dx dP] = gradient(P(:,:,end));
    dPV = sum(dP.*dP,2).^.5;
    fidx = find(dPV > THRESH);
    toTrack(fidx) = 0;            
    for f = 1:numel(fidx)
        if fidx(f) ~= 1
            if fidx(f) ~= size(P,1)
                dP(fidx(f),:,:) = mean(dP(fidx(f)-1:fidx(f)+1,:,:),1);
                plot(P(p,1,e+1),P(p,2,e+1),'r*')
            end
        end
    end
    %}
            
end
%%
%[c r V] = impixel(to);
close all
c = (1:1:size(I,2))';
%c = linspace(200,size(I,2),2000)';
r = (size(I,1)-100*ones(1,size(c,1)))';
P = [c,r];
THRESH = 10;
toTrack = ones(1,size(P,1));
imshow(to);hold on
for e = 1:200
    for p = 1:size(P,1)
        if toTrack(p)
            d1 = ba_interp2(W1,P(p,1,e),P(p,2,e));
            d2 = ba_interp2(W2,P(p,1,e),P(p,2,e));
            d = [d1 d2];
            d = d / norm(d);
            d = d * .5;

            if e ~= 1
                dV = squeeze(P(p,:,e) - P(p,:,e-1));
                if d*dV'< 0
                    d = -d;                
                end
            end

            P(p,1,e+1) = P(p,1,e) + d(1);
            P(p,2,e+1) = P(p,2,e) + d(2);
        end
        p
        e
        
    end
        
        
        
     %{   
    [dx dP] = gradient(P(:,:,end));
    dPV = sum(dP.*dP,2).^.5;
    fidx = find(dPV > THRESH);
    toTrack(fidx) = 0;            
    for f = 1:numel(fidx)
        if fidx(f) ~= 1
            if fidx(f) ~= size(P,1)
                dP(fidx(f),:,:) = mean(dP(fidx(f)-1:fidx(f)+1,:,:),1);
                plot(P(p,1,e+1),P(p,2,e+1),'r*')
            end
        end
    end
    %}
            
    e
    size(P,1)
end


for p = 1:size(P,1)
    if toTrack(p)
        plot(squeeze(P(p,1,:)),squeeze(P(p,2,:)));
    end
end

drawnow

%% find files and generate rand number
pth = '/mnt/scratch3/users/nmiller/phytoMorph/outPort/'
dataFiles = gdig(pth,{},{'mat'},1);
numel(dataFiles);
%%
parfor e = 1:100
    d = load(dataFiles{e})
    d = d.d;
    rootN = 1;
    I = imread(d{rootN}{1}.image.fileName);    
    %imshow(I,[])
    %drawnow
    try
    S{e} = I;
    catch
    end
    e
    
end
%%
Tr.learn();
%%
Tr.viewManifold();
%% patch

h = figure;
patchS.view(h,[],[]);

%% look at curves

close all

Tr.viewCoandDo();
close all




