FilePath = '/mnt/snapper/nate/mirror_images/DRB/Gravatron_PS/';
FilePath = '/mnt/snapper/nate/drB/FL_4_with_box_QRcode_v1';
FileList = {};
FileExt = {'jpg'};
FileList = sdig(FilePath,FileList,FileExt,1);
%%
for s = 1:numel(FileList)
    tmp = imread(FileList{s});
    I = zeros([size(tmp) numel(FileList{s})]);
    for e = 1:numel(FileList{s})
        I(:,:,:,e) = imread(FileList{s}{e});
    end
    %% make gray
    for e = 1:size(I,4)
        G(:,:,e) = double(rgb2gray(I(:,:,:,e)));
        fprintf(['Done with gray for:' num2str(e) ':' num2str(size(I,4)) '\n']);
    end
    %% whole thing
    [sidx] = findRotateFrame(G);
    [bk dM] = getBKandDM(I(:,:,:,1:sidx),20);
    sd = std(dM,1,3);
    %{
    %% background
    BK = mean(double(G(:,:,1:10)),3);
    imshow(BK,[])
    %% background movie
    close all
    d = bsxfun(@minus,double(G(:,:,1:sidx)),BK);
    sd = std(d,1,3);
    %}
    %% analysis of sd
    sd1 = mean(sd,2);
    [~,tidx] = max(sd1);
    sd2 = mean(sd(tidx:tidx+200,:),1);
    sd2 = bindVec(sd2);
    cidx = find(imdilate(sd2,strel('disk',100)) == sd2);
    bsig = sd2 > graythresh(sd2);
    cidx(~bsig(cidx)) = [];
    %% show line for box
    WINBOX(1) = [75];
    close all
    tidx = tidx - WINBOX(1)
    imshow(I(:,:,:,sidx),[]);
    hold on
    plot(1:size(I,2),tidx*ones(1,size(I,2)),'r')
    for e = 1:numel(cidx)
        plot(cidx(e),tidx,'r*');
        tBOX = [cidx(e)-100 tidx 200 1000];
        tmpI = imcrop(I(:,:,:,1),tBOX);
        imshow(tmpI,[]);
        drawnow
        BOX{e} = tBOX;
        cI{e} = cropMovie(I(:,:,:,1:sidx),tBOX);
    end
    %% process seedling stack
    [SHOOT,ROOT] = processSeedlingStacks(cI);
    %% view on raw data
    close all
    writerObj = VideoWriter('/mnt/snapper/nate/mirror_images/DRB/mov1.avi') ;
    writerObj.open();
    for t = 1:sidx
        imshow(I(:,:,:,t),[]);
        hold on
        for s = 1:numel(SHOOT)
            plot(SHOOT(s).X(t,:)+BOX{s}(1),SHOOT(s).Y(t,:)+BOX{s}(2),'g')
            plot(ROOT(s).X(t,:)+BOX{s}(1),ROOT(s).Y(t,:)+BOX{s}(2),'r')
        end
        drawnow
        frame = getframe;
        writeVideo(writerObj,frame);
        hold off
    end
    writerObj.close();
    
end
%% make gray
for e = 1:size(I,4)
    G(:,:,e) = double(rgb2gray(I(:,:,:,e)));
    fprintf(['Done with gray for:' num2str(e) ':' num2str(size(I,4)) '\n']);
end
%% whole thing
[sidx] = findRotateFrame(G);
[bk dM] = getBKandDM(I(:,:,:,1:sidx),20);
sd = std(dM,1,3);
%{
%% background
BK = mean(double(G(:,:,1:10)),3);
imshow(BK,[])
%% background movie
close all
d = bsxfun(@minus,double(G(:,:,1:sidx)),BK);
sd = std(d,1,3);
%}
%% analysis of sd
sd1 = mean(sd,2);
[~,tidx] = max(sd1);
sd2 = mean(sd(tidx:tidx+200,:),1);
sd2 = bindVec(sd2);
cidx = find(imdilate(sd2,strel('disk',100)) == sd2);
bsig = sd2 > graythresh(sd2);
cidx(~bsig(cidx)) = [];
%% show line for box
WINBOX(1) = [75];
close all
tidx = tidx - WINBOX(1)
imshow(I(:,:,:,sidx),[]);
hold on
plot(1:size(I,2),tidx*ones(1,size(I,2)),'r')
for e = 1:numel(cidx)
    plot(cidx(e),tidx,'r*');
    tBOX = [cidx(e)-100 tidx 200 1000];
    tmpI = imcrop(I(:,:,:,1),tBOX);
    imshow(tmpI,[]);
    drawnow
    BOX{e} = tBOX;
    cI{e} = cropMovie(I(:,:,:,1:sidx),tBOX);
end
%% process seedling stack
[SHOOT,ROOT] = processSeedlingStacks(cI);
%% view on raw data
close all
writerObj = VideoWriter('/mnt/snapper/nate/mirror_images/DRB/mov1.avi') ;
writerObj.open();
for t = 1:sidx
    imshow(I(:,:,:,t),[]);
    hold on
    for s = 1:numel(SHOOT)
        plot(SHOOT(s).X(t,:)+BOX{s}(1),SHOOT(s).Y(t,:)+BOX{s}(2),'g')
        plot(ROOT(s).X(t,:)+BOX{s}(1),ROOT(s).Y(t,:)+BOX{s}(2),'r')
    end
    drawnow
    frame = getframe;
    writeVideo(writerObj,frame);
    hold off
end
writerObj.close();
%%
sm = [];
for r = 1:size(I,4)
    sm(:,:,:,r) = imresize(I(:,:,:,r),.10);
end
%%
sz = size(sm);
bk = double(reshape(sm,[prod(sz(1:3)) sz(4)]));
[S C U E L ERR LAM] = PCA_FIT_FULL_T(bk(:,1:10),3);
