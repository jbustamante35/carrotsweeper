%% 1) scan for mat
FilePath = '/mnt/tetra/nate/chips2/';
FileList = {};
FileExt = {'mat'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% 2) load
toRM = [];
STACK = zeros(1500,3,60,numel(FileList));
GRIDY = zeros(1500,2,60,numel(FileList));
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
%% backup
STACK_BK = STACK;
GRIDY_BK = GRIDY;
%% restore
STACK = STACK_BK;
GRIDY = GRIDY_BK;
%% align
RM = [];
for tr = 1:size(GRIDY,4)    
    for ls = 1:size(GRIDY,3)
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
    fprintf(['done aligning:' num2str(tr) ':' num2str(size(STACK,4)) '\n']);
end
%% remove from align
GRIDY(:,:,:,RM) = [];
STACK(:,:,:,RM) = [];
IG(RM) = [];
%% plot size shape
figure;
plot(squeeze(mean(GRIDY(:,2,1,:),4)),squeeze(mean(GRIDY(:,1,1,:),4)),'r')
hold on
plot(squeeze(mean(GRIDY(1,2,1,:),4)),squeeze(mean(GRIDY(1,1,1,:),4)),'r*')
%% make more
STACK = cat(4,STACK,flipdim(STACK,1));
%% make sym
STACK = .5*(STACK + flipdim(STACK,1));
fprintf(['done SYM:' num2str(tr) ':' num2str(size(STACK,4)) '\n']);
%% rgb to Lab
STACK = ipermute(rgb2lab(permute(STACK,[1 3 2 4])),[1 3 2 4]);
fprintf(['done rgb->Lab:' num2str(tr) ':' num2str(size(STACK,4)) '\n']);
%% try ring image
close all
Rmin = 50;
EXTRA = 0;
Rmax = Rmin + size(STACK,3) + EXTRA;
[g1 g2] = ndgrid(linspace(Rmax,Rmin,size(STACK,3)),linspace(-pi,pi,size(STACK,1)));
G1 = g1.*cos(g2-pi/2);
G2 = g1.*sin(g2-pi/2);
h1 = figure;
h2 = figure;
for toShow = 1:10
    [Z] = makeRingImage(STACK(:,:,:,toShow));
    
    figure(h1);
    imshow(Z(:,:,1),[])
    drawnow
    
    figure(h2);
    imgTmp = IG{toShow};
    imgTmp = rgb2lab(imgTmp);
    imshow(imgTmp(:,:,2),[])
    
    waitforbuttonpress
end

%% non fft
sz = size(STACK);
STACK = reshape(STACK,[prod(sz(1:3)) sz(4)]);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL_T(STACK,4);
%% sweep for non fft - with new ring
close all
U = squeeze(reshape(sU,[sz(1:(end-1))]));
imshow(U,[]);
NS = 10;
US = [];
for nc = 2%1:size(sC,1)
    MU = mean(sC,2);
    L = linspace(min(sC(nc,:)),max(sC(nc,:)),10);
    for l = 1:numel(L)
        tmpC = MU;
        tmpC(nc) = L(l);
        tmp = PCA_BKPROJ_T(tmpC,sE,sU);
        U = (reshape(tmp,[sz(1:(end-1))]));
        U = imresize(squeeze(U),[size(U,1) 3*size(U,3)]);
        U = reshape(U,[size(U,1) 1 size(U,2)]);
        US(:,:,l) = makeRingImage(U);        
    end
    
    for loop = 1:3
        for l = 1:numel(L)
            imshow(US(:,:,l),[])
            title([num2str(nc) '--' num2str(l)]);
            drawnow
            pause(.3)
        end
    end
end
%% isolate channal - 
STACK = STACK(:,2,:,:);
fprintf(['done isolate channel:' num2str(tr) ':' num2str(size(STACK,4)) '\n']);
%% smooth
for e3 = 1:size(STACK,4)
    for e1 = 1:size(STACK,2)
        for e2 = 1:size(STACK,3)        
            STACK(:,e1,e2,e3) = imfilter(STACK(:,e1,e2,e3),fspecial('average',[51 1]),'circular');
        end
    end
    fprintf(['done smoothing:' num2str(e3) ':' num2str(size(STACK,4)) '\n']);
end
STACK = permute(STACK,[3 1 2 4]);
for e4 = 1:size(STACK,4)
    for e3 = 1:size(STACK,3)
        STACK(:,:,e3,e4) = imfilter(STACK(:,:,e3,e4),fspecial('average',[3 1]),'replicate');
    end
    fprintf(['done smoothing:' num2str(e4) ':' num2str(size(STACK,4)) '\n']);
end
STACK = ipermute(STACK,[3 1 2 4]);
%% reshape
STACK = reshape(STACK,sz);
%% diff - along rad
[d1 d2 d3 d4] = gradient(STACK);
STACK = (d2.^2 + d3.^2).^.5;
%STACK = d2;
%STACK = (abs(diff(STACK,1,1)).^2 + abs(diff(STACK,1,3)).^2).^.5;
%% diff - along theta
[d1 d2 d3 d4] = gradient(cat(1,flipdim(STACK(1:10,:,:,:),1),STACK,flipdim(STACK((end-9:end),:,:,:),1)));
%STACK = (d2.^2 + d3.^2).^.5;
STACK = d2;
STACK = STACK(11:(end-10),:,:,:);
%STACK = (abs(diff(STACK,1,1)).^2 + abs(diff(STACK,1,3)).^2).^.5;
%% post process backup
HSTACK = STACK;
%% post process restore
STACK = HSTACK;
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
SKIP = 700;
STOP = 6000;
pt = [];
tightSTRING = 'COM';
for r = 1:size(sC,1)
    
MQ = [];
MQ2 = [];
MQ3 = [];

    
    close all
    h1 = figure;
    h2 = figure;
    toSave = 0;
    mon = 1;
    oPath = '/mnt/snapper/nate/mirror_images/forAnalysis/return/FFT_MORE/';
    GG = 1:4;
  
    
    
    GG = setdiff(GG,r);


    pth = zeros(size(sC,1),10);
    
    pth(r,:) = linspace(min(sC(r,:)),max(sC(r,:)),10);
    for str = 1:3
        for p = 1:size(pth,2)
            d = bsxfun(@minus,sC,pth(:,p));
            [d sidx] = sort(sum(d.*d,1));
            for l = 1:numel(pt)
                sidx(sidx==pt(l)) = [];
            end
            %sidx = setdiff(sidx,pt);
            pt(str,p) = sidx(1);
            tmpsC(:,str,p) = sC(:,sidx(1));
        end
    end
    
    
    for str = 1:3
    
      
        Q = [];
        Q2 = [];
        Q3 = [];
    
        sidx = pt(str,:);
        
        
        %{
        
        [J sidx] = sort(sC(r,:));
        
        %J = abs(sC) < 7000;
        J = abs(sC) < 500/1;
        DT = sum(sC(GG,sidx).*sC(GG,sidx),1);
        
        
        [DT,fil] = sort(DT);
        
        %fil = find(all(J(GG,sidx),1));
        
        %sidx = sidx(fil(1:100));
        
        %sidx = pt;
        
        %[J sidx] = sort(Out1(r,:));
        %ssS = sS(:,sidx);
        %}
        %[J sidx] = sort(sC(r,:));
        nsz = [500 500];
      
        tP = [oPath num2str(r) filesep];
        mkdir(tP)
        
        %for e = str:SKIP:STOP+str
        for e = 1:numel(sidx)%1:round(numel(sidx)/10):numel(sidx)
            %{
            imshow(IGF{sidx(e)} ,[]);
            drawnow
            waitforbuttonpress
            %}
            tmpIMG = IG{sidx(e)};
            tmpIMG_LAB = rgb2lab(tmpIMG);
            %tmpIMG_LAB = tmpIMG;
            figure(h1);
            fileName = [tP num2str(e) '.jpg'];
            if toSave
                imwrite(tmpIMG,fileName);
            end
            imshow(tmpIMG_LAB(:,:,2),[-5 20]);
            title([num2str(r) '--' num2str(e)]);
            drawnow
            %M(r,e,:) = size(IGF{sidx(e)});
            tmpC = tmpsC(:,str,e);
            tmpC(GG) = 0;
            tmpOver = PCA_BKPROJ_T(tmpC,sE,sU);
            tmpOver= (reshape(tmpOver,[sz(1:(end-1))]));
            tsz = [size(GRIDY,1) 1 size(GRIDY,3)];
            g1 = reshape(imresize(squeeze(GRIDY(:,1,:,sidx(e))),[tsz(1)*2 tsz(3)*2]),[tsz(1)*2 1 tsz(3)*2]);
            g2 = reshape(imresize(squeeze(GRIDY(:,2,:,sidx(e))),[tsz(1)*2 tsz(3)*2]),[tsz(1)*2 1 tsz(3)*2]);
            tmpOver = reshape(imresize(squeeze(tmpOver(:,1,:)),[tsz(1)*2 tsz(3)*2]),[tsz(1)*2 1 tsz(3)*2]);
            
            tGRID = cat(2,g1,g2);            
            tmpOver = makeRingImage_org(tmpIMG,tmpOver,tGRID);


            
            %{
            figure(h2)
            plot(sS(:,sidx(e)));
            axis([0 size(STACK,1) 0 max(sS(:))]);
            drawnow
            %}
            if mon 
                Q = cat(2,Q,imresize(tmpIMG_LAB,nsz));
                Q2 = cat(2,Q2,imresize(tmpOver,nsz));
                Q3 = cat(2,Q3,imresize(tmpIMG,nsz));
            end
        end
        %figure;imshow(Q,[]);
        %waitforbuttonpress
        
        
        MQ = cat(1,MQ,Q);
        MQ2 = cat(1,MQ2,Q2);
        MQ3 = cat(1,MQ3,Q3);
        
    end
    imwrite(bindVec(MQ(:,:,2)),['/mnt/tetra/' tightSTRING '--' num2str(r) '--Lab.tif']);
    imwrite(bindVec(MQ2),['/mnt/tetra/' tightSTRING '--' num2str(r) '--PCA.tif']);    
    imwrite(MQ3,['/mnt/tetra/' tightSTRING '--' num2str(r) '--RGB.tif']);
end
%%
close all
imshow(MQ(:,:,2),[])
figure;
imshow(MQ2,[])
figure;
imshow(MQ3,[])
imwrite(MQ,'/mnt/tetra/kurt.tif')
%%
close all
imshow(MQ(:,:,:),[])
imwrite(MQ,'/mnt/tetra/kurt.tif')
