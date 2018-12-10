%% 1) scan for mat
FilePath = '/mnt/tetra/nate/chips/';
FileList = {};
FileExt = {'mat'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% 2) load
toRM = [];
STACK = zeros(1000,3,25,numel(FileList));
GRIDY = zeros(1000,2,25,numel(FileList));
IG = cell(1,numel(FileList));
for e = 1:numel(FileList)
    try
        tmp = load(FileList{e},'sam','img','grid');
        STACK(:,:,:,e) = tmp.sam;
        GRIDY(:,:,:,e) = tmp.grid;
        IG{e} = tmp.img;
        
    catch
        toRM = [toRM e];
    end
    fprintf(['done loading:' num2str(e) ':' num2str(numel(FileList)) '\n']);
end
STACK(:,:,:,toRM) = [];
GRIDY(:,:,:,toRM) = [];
IG(toRM) = [];
%%
STACK = [];
GRIDY = [];
%IGF = {};
for e = 1:numel(FUN)
    if ~isempty(FUN)        
        STACK = cat(4,STACK,FUN{e});
        GRIDY = cat(4,GRIDY,GRID{e});
        %{
        for l = 1:numel(IG{e})
            IGF{end+1} = rgb2lab(IG{e}{l});
        end
        %}
    end
end
%% backup
STACK_BK = STACK;
GRIDY_BK = GRIDY;
%% restore
STACK = STACK_BK;
GRIDY = GRIDY_BK;
%% stack images
cnt =1;
IGF = cell(size(STACK,4),1);
for e = 1:numel(FUN)
    if ~isempty(FUN) 
        for l = 1:numel(IG{e})
            IGF{cnt} = rgb2lab(IG{e}{l});
            cnt = cnt + 1;
        end
    end
end
%% align
RM = [];
for ls = 1:size(GRIDY,3)
    for tr = 1:size(GRIDY,4)
        try
                fidx = find(GRIDY(:,2,ls,tr) > mean(GRIDY(:,2,ls,tr)));
                %{
                plot(GRID{e}(fidx,2,ls,tr),GRID{e}(fidx,1,ls,tr))
                drawnow
                waitforbuttonpress
                %}
                STACK(:,:,ls,tr) = circshift(STACK(:,:,ls,tr),-fidx(end),1);
                GRIDY(:,:,ls,tr) = circshift(GRIDY(:,:,ls,tr),-fidx(end),1);
        catch
            RM = [RM tr];
        end
    end
end
%% remove
GRIDY(:,:,:,RM) = [];
STACK(:,:,:,RM) = [];
%% plot size shape
figure;
plot(squeeze(mean(GRIDY(:,2,1,:),4)),squeeze(mean(GRIDY(:,1,1,:),4)),'r')
hold on
plot(squeeze(mean(GRIDY(1,2,1,:),4)),squeeze(mean(GRIDY(1,1,1,:),4)),'r*')
%% make more
STACK = cat(4,STACK,flipdim(STACK,1));
%% rgb to Lab
STACK = ipermute(rgb2lab(permute(STACK,[1 3 2 4])),[1 3 2 4]);
%% non fft
sz = size(STACK);
STACK = reshape(STACK,[prod(sz(1:3)) sz(4)]);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL_T(STACK,4);
%% sweep for non fft
U = squeeze(reshape(sU,[sz(1:(end-1))]));
imshow(U,[]);
for nc = 1:size(sC,1)
    MU = mean(sC,2);
    L = linspace(min(sC(nc,:)),max(sC(nc,:)),5);
    for l = 1:numel(L)
        tmpC = MU;
        tmpC(nc) = L(l);
        tmp = PCA_BKPROJ_T(tmpC,sE,sU);
        U = squeeze(reshape(tmp,[sz(1:(end-1))]));
        imshow(imresize(U,[500 500]),[])
        title([num2str(nc) '--' num2str(l)]);
        drawnow
        waitforbuttonpress
    end
end
%% isolate channal - 
STACK = STACK(:,2,:,:);
%% smooth
for e1 = 1:size(STACK,2)
    for e2 = 1:size(STACK,3)
        for e3 = 1:size(STACK,4)
            STACK(:,e1,e2,e3) = imfilter(STACK(:,e1,e2,e3),fspecial('average',[51 1]),'circular');
        end
    end
end
%% diff
STACK = (abs(diff(STACK,1,1)).^2 + abs(diff(STACK,1,3)).^2).^.5;
%% mean along chip
STACK = squeeze(mean(STACK,3));
%% fft
ND = size(STACK,1)/2;
ND = 100;
%STACK = bsxfun(@minus,mean(STACK,1),STACK);
%STACK = bsxfun(@minus,mean(STACK,3),STACK);
%STACK = bsxfun(@times,STACK,sum(STACK.*STACK,1).^-.5);
STACK = abs(fft(STACK,[],1));
STACK = STACK(1:ND,:,:,:);
CP = STACK;
%FY = squeeze(mean(mean(STACK,3),4));
sz = size(STACK);
STACK = reshape(STACK,[prod(sz(1:(end-1))) sz(end)]);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL_T(STACK,4);
%% sweep along fft
%%
MQ = [];
%%
MQ = [];

SKIP = 60;
STOP = 600;

for str = 1:4
    close all
    h1 = figure;
    h2 = figure;
    toSave = 0;
    mon = 1;
    oPath = '/mnt/snapper/nate/mirror_images/forAnalysis/return/FFT_MORE/';
    for r = 1%:size(sC,1)
        [J sidx] = sort(sC(r,:));
        %[J sidx] = sort(Out1(r,:));
        ssS = sS(:,sidx);
        nsz = [500 500];
        Q = [];
        tP = [oPath num2str(r) filesep];
        mkdir(tP)
        
        for e = str:SKIP:STOP+str
            %{
            imshow(IGF{sidx(e)} ,[]);
            drawnow
            waitforbuttonpress
            %}
            tmpIMG = IG{sidx(e)};
            tmpIMG_LAB = rgb2lab(tmpIMG);
            figure(h1);
            fileName = [tP num2str(e) '.jpg'];
            if toSave
                imwrite(tmpIMG,fileName);
            end
            imshow(tmpIMG_LAB(:,:,2),[-5 20]);
            title([num2str(r) '--' num2str(e)]);
            drawnow
            %M(r,e,:) = size(IGF{sidx(e)});



            %{
            figure(h2)
            plot(sS(:,sidx(e)));
            axis([0 size(STACK,1) 0 max(sS(:))]);
            drawnow
            %}
            if mon 
                Q = cat(2,Q,imresize(tmpIMG_LAB,nsz));
            end
        end
        %figure;imshow(Q,[]);
        %waitforbuttonpress
    end
    MQ = cat(1,MQ,Q);
end
%close all
imshow(MQ(:,:,2),[])
%%
close all
imshow(MQ(:,:,2),[])
