FilePath = '/mnt/spaldingdata/nate/mirror_images/04202015 Arabidopsis/04202015 mock/';
FilePath = '/mnt/spaldingdata/nate/mirror_images/04272015 tomato/04272015 tomato mock/';
FilePath ='/mnt/spaldingdata/nate/mirror_images/04242015 Arabidopsis/iaa/';
FilePath ='/mnt/spaldingdata/nate/mirror_images/04242015 Arabidopsis/mock/';
miniP = '05052015 dim light/05052015 Arabidopsis dim light iaa/';
miniP = '05052015 dim light/05052015 Arabidopsis dim light mock/';
miniP = '05052015 etiolated/05052015 etiolated mock/';
miniP = '05052015 tomato/05052015 tomato mock/';
miniP = '05052015 tomato/05052015 tomato iaa/';
FilePath = '/mnt/spaldingdata/nate/mirror_images/forBill/';
FilePath = [FilePath miniP];
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
FileList(1:2) = [];
%% sort the images
n = [];
for e = 1:numel(FileList)
    [p,tn,ex] = fileparts(FileList{e});
    n(e) = str2num(tn);    
end
[~,sidx] = sort(n);
FileList = FileList(sidx);
%% read the images
for e = 1:numel(FileList)
    I(:,:,e) = imread(FileList{e});    
    e
end
%% view the movie
close all
for e = 1:numel(FileList)
    imshow(I(:,:,e),[]);
    drawnow
end
%% binary the images
B = [];
S = [];
parfor e = 1:numel(FileList)
    tmp = double(I(:,:,e))/255;
    tmp = imfilter(tmp,fspecial('gaussian',[21 21],9),'replicate');
    level = graythresh(tmp);
    B(:,:,e) = double(tmp) < level;
    B(:,:,e) = bwareaopen(B(:,:,e),1000);
    B(:,:,e) = imfill(B(:,:,e),'holes');
    B(:,:,e) = imclearborder(B(:,:,e));
    S(:,:,e) = bwmorph(B(:,:,e),'skel',inf);
    tmp = S(:,:,e);
    for k = 1:20
        ep = bwmorph(tmp,'endpoints');
        tmp(find(ep)) = 0;
    end
    S(:,:,e) = tmp;
    e
end
%% get contours and midlines
%dB = cell(numel(FileList));
P = cell(numel(FileList));
pidx = cell(numel(FileList));
loc = cell(numel(FileList));

parfor e = 1:numel(FileList)    
    dB{e} = bwboundaries(B(:,:,e));
    for c = 1:numel(dB{e})
        W = poly2mask(dB{e}{c}(:,2), dB{e}{c}(:,1), size(B,1),size(B,2));
        W = S(:,:,e).*W;
        [iy ix] = find(W);
        [out] = cwtK_closed_imfilter(dB{e}{c},{[51]});
        sig = [out.K;out.K;out.K];
        peaks = imdilate(-sig,strel('disk',50,0)) == -sig;
        str = size(out.K,1)+1;
        stp = str + size(out.K) -1;
        peaks = peaks(str:stp);
        pidx{e}{c} = find(peaks);
        values = out.K(pidx{e}{c});
        [~,sidx] = sort(values);
        pidx{e}{c} = pidx{e}{c}(sidx(1:2));
        loc{e}{c} = dB{e}{c}(pidx{e}{c},:);         
        DP = [fliplr(loc{e}{c})' [ix iy]'];
        T = Radjacency(DP,100);
        [path , pathcost]  = dijkstra(T , 1 , 2);
        P{e}{c} = DP(:,path);
        fprintf(['Done with analysis: ' num2str(c) ' and image ' num2str(e) '\n']);
    end
end
%% show binary movie
close all
for e = 1:numel(FileList)    
    imshow(B(:,:,e),[]);
    drawnow
end
%% show midlines and contours 
figure;
%close all
Z = zeros(size(B,1),size(B,2));
for e = 1:numel(FileList)    
    RGB = cat(3,S(:,:,e),Z,B(:,:,e));        
    imshow(double(RGB),[]);
    hold on
    for c = 1:numel(dB{e})
        plot(dB{e}{c}(:,2),dB{e}{c}(:,1),'g')
        plot(loc{e}{c}(:,2),loc{e}{c}(:,1),'r*')
        plot(P{e}{c}(1,:),P{e}{c}(2,:),'y')
    end
    hold off
    title(num2str(e))
    drawnow
end
%% measure
L = [];
WID = [];
for img = 1:numel(FileList)
    tmp = bwdist(~B(:,:,img));
    for con = 1:numel(dB{e})
        
        IDX = sub2ind(size(tmp),P{img}{con}(2,:),P{img}{con}(1,:));
        values = tmp(IDX);
        WID(img,con) = mean(values);
        dL = diff(P{img}{con},1,2);
        dL = sum(dL.*dL,1).^.5;
        L(img,con) = sum(dL);
        LEG{con} = num2str(con);
    end
    img
end
%% plot measurements
close all
plot(L)
legend(LEG)
figure
plot(WID)
legend(LEG)
figure;
e = 1;
RGB = cat(3,S(:,:,e),Z,B(:,:,e));        
imshow(double(RGB),[]);
hold on
for c = 1:numel(dB{e})
    plot(dB{e}{c}(:,2),dB{e}{c}(:,1),'g')
    plot(loc{e}{c}(:,2),loc{e}{c}(:,1),'r*')
    plot(P{e}{c}(1,:),P{e}{c}(2,:),'y')
    uP = mean(dB{e}{c},1);
    text(uP(2),uP(1),num2str(c),'BackgroundColor',[1 1 1]);
end
hold off
drawnow
%%


