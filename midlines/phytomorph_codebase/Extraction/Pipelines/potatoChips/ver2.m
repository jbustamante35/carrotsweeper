%% 1) scan for images
FilePath = '/mnt/spaldingdata/nate/mirror_images/PotatoChips/';
FilePath = '/mnt/snapper/nate/mirror_images/forAnalysis/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%%
for e = 1:3:numel(FileList)
    info = imfinfo(FileList{e});
    EL = {round([1 10 info.Height]) round([1 10 info.Width])}
    tmp = imread(FileList{e},'PixelRegion',EL);
    imshow(tmp,[])
    waitforbuttonpress
end
%% 2) get 2015 images only
kidx = [];
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'2015'))
        kidx= [kidx e];
    end
end
FileList = FileList(kidx);
%%
for e = 1:N
    imf = imfinfo(FileList{e});
    imf.Width
    imf.Height
end
%% X) view some images from 1 -> N
IS = [];
N = 10;
for e = 1:N
    I = imread(FileList{e});
    %imshow(I,[])
    
    G = rgb2gray(I);
    imshow(I,[]);
    title(size(G))
    drawnow
    try
        IS = cat(3,IS,G);
    catch ME
        ME
    end
    size(IS)
    e
end
%% NOT NEEDED make mask
I = imread(FileList{1});
Z = zeros(size(I,1),size(I,2));
Z(1:1300,:) = 1;
%% NOT NEEDE try gradient
for e = 1:numel(FileList)
    I = imread(FileList{e});
    G = double(rgb2gray(I));
    [g1 g2] = gradient(G);
    dG = (g1.^2 + g2.^2).^.5;
end
%% 4) try threshold again - as of May 21 2015 this seems to work
% gets qr msg and masks l-a-b
close all
disp = 1;
msg = {};
for e = 1:10%numel(FileList)
    try
        tic
        %if isempty(msg{e})
            % get the chip mask
            [I B] = getChipMask_ver2(FileList{e},5000,1000);
            [qr msg{e}] = getQR(I,B);
            %{

            [QR] = imcrop(I);
            QR = rgb2gray(QR);
            QR = flipdim(QR,1);
            msg = decode_qr(QR);
            %}
            % convert the RGB to LAB
            SAM = [];
            fidx = find(B);
            for k = 1:size(I,3)
                tmp = I(:,:,k);
                SAM(:,k) = tmp(fidx);
            end
            % apply conversion
            C = makecform('srgb2lab');
            LAB = applycform(SAM/255,C);
            %{
            lab = zeros(size(I));
            fidx = find(B);
            for k = 1:size(I,3)
                tmp = zeros(size(I,1),size(I,2));
                tmp(fidx) = LAB(:,k);
                lab(:,:,k) = tmp;
            end
            %}
            % get the average LAB
            uLAB(e,:) = mean(LAB,1);


           
            if disp
                 % get the boundaries for the mask
                dB = bwboundaries(B);
                % show the boundaries and the image
                close all
                imshow(I,[]);
                hold on    
                for i = 1:numel(dB)
                    plot(dB{i}(:,2),dB{i}(:,1),'r')
                end


                hold off
                title(num2str(e));
                drawnow
            end
        %end
        toc
        e
    catch ME
        msg{e} = [];
        fprintf(['error in:' num2str(e) '\n']);
    end
end
%%
[J1,J2,caliData] = xlsread('/home/nate/Downloads/2014-15_chips_hunter.xlsx',1);
[J1,J2,mLAB] = xlsread('/home/nate/Downloads/2014-15_chips_hunter.xlsx',24);
for e = 2:size(caliData,1)
   
end
%% get qr msg only and after
msg = {};
for e = 1:numel(FileList)
    % get the chip mask
    I = imread(FileList{e});
    [qr msg{e}] = getQR(I);
end
%% focus on stem end defect with data mining and parameterization of color space
% note the first image is as a different resolution
close all
tcH = cell(1,100);
disp = 0;
% make n hood to sample the stem-end which is "always" up
parfor e = 1:numel(FileList)
    
    %imgInfo = imfinfo(FileList{e});
    %if imgInfo.Width < 4000
        e
        % read the image
        I = imread(FileList{e});
        % get the chip mask
        [I B] = getChipMask_ver0(FileList{e},5000,1000);
        B = logical(B);
        fprintf(['Done with mask \n']);
        % make inner ship mask and then subtract from chip mask to reveal the
        % outter rim
        h = strel('disk',51,0);
        iB = imerode(B,h);        
        fprintf(['Done with erode ring \n']);
        % get the boundaries for the mask
        idB = bwboundaries(iB);
        % create the ring mask
        rM = B - iB;
        % sample the ring mask
        %HI = sampleRing(rM,I);
        fprintf(['Done with sample ring \n']);
        % make the total count
        %tcH = cat(3,tcH,HI);
        tcH{e} = HI;

        R = regionprops(B,'PixelIdxList','BoundingBox');
        % show the user the chips in order and ask a question


        
        
        % get the boundaries for the mask
        dB = bwboundaries(B);
        
        [FUN{e} GRID{e} IG{e}] = ringSample(I,B,dB,R);
        
        h1 = figure;
        h2 = figure;
        
        
        if disp
            figure(h1);
            % show the boundaries and the image            
            imshow(I,[]);
            hold on    
            for i = 1:numel(dB)
                plot(dB{i}(:,2),dB{i}(:,1),'r')
            end
            for i = 1:numel(dB)
                plot(idB{i}(:,2),idB{i}(:,1),'g')
            end
            
            for i = 1:numel(dB)
                plot(mean(idB{i}(:,2)),mean(idB{i}(:,1)),'b*')
            end
            
            
            
            
            
            title(num2str(e));
            drawnow
            
            
            ORD = [];
            
            for i = 1:numel(dB)
                for k = 1:3
                    ORD(i,k) = mean(-log(HI(1:10,k,i)),1);
                end
            end
            
            CL = [-150 0 150];
            CLL = {'r' 'g' 'b'};
            for k = 1:3
                [J sidx] = sort(ORD(:,k));
                for i = 1:numel(sidx)
                    tmp = mean(idB{sidx(i)},1);
                    text(tmp(2)-CL(k),tmp(1),num2str(i),'BackGroundColor',CLL{k});
                end
            end
            
            %{
            for i = 1:numel(dB)
                figure(h1)
                plot(idB{i}(:,2),idB{i}(:,1),'b')
                figure(h2);
                plot(-log(HI(:,:,i)))
                waitforbuttonpress
            end
            %}
            
        end
        
        
    %end
    e
   
end
%%
parfor e = 1:numel(FileList)
    er(e) = pp_ver0(FileList{e});
end

%% stack the tcH
H = [];
NM = [];
for e = 1:numel(tcH)
    if size(tcH{e},1) ~= 0
        H = cat(3,H,tcH{e});
        NM = [NM;e*ones(size(tcH{e},3),1)];
    end
end
%% reshape
sz = size(H);
H = reshape(H,[prod(sz(1:2)) sz(3)]);
%% show the 3D plot for color space
[S C U E L ERR LAM] = PCA_FIT_FULL_T(H,3);
plot3(C(1,:),C(2,:),C(3,:),'.');
%% peform 5 kmeans groups on the data - CAREFUL
close all
kidx = kmeans(C',4);
UQ = unique(kidx);
CL = {'r.' 'g.' 'b.' 'c.'};
for u = 1:numel(UQ)
    fidx = kidx==u;
    plot3(C(1,fidx),C(2,fidx),C(3,fidx),CL{u});
    hold all
end
%% permute the data - careful
kkidx = kidx;
kp = [2 4 1 3];
for e = 1:4
    fidx = find(kkidx==kp(e));
    kidx(fidx) = e;
end
%% display the permutations
close all
UQ = unique(kidx);
CL = {'r.' 'g.' 'b.' 'c.'};
for u = 1:numel(UQ)
    fidx = kidx==u;
    plot3(C(1,fidx),C(2,fidx),C(3,fidx),CL{u});
    hold all
end
%% smooth parameter?
UQ = unique(kidx);
for u = 1:numel(UQ)
    fidx = kidx==u;
    uC(:,u) = mean(C(:,fidx),2)
end
dU = diff(uC,1,2);
dU = cumsum([0 sum(dU.*dU,1).^.5]);
dU = dU / max(dU);
spl = spap2(1,3,dU,uC);
Xi = linspace(0,max(dU),1000);
SM = fnval(spl,Xi);
plot3(SM(1,:),SM(2,:),SM(3,:),'k');
%% order the points along the curve
close all
for e = 1:size(C,2)
    delta = bsxfun(@minus,SM,C(:,e));
    [J,loc(e)]= min(sum(delta.*delta,1));
end
loc = Xi(loc);
[J,sidx] = sort(loc);
dist = .1;
eE = [];
for e = 1:numel(Xi)
    fidx = find(abs(loc-Xi(e)) <= dist);
    plot3(SM(1,:),SM(2,:),SM(3,:),'g')    
    hold on
    plot3(SM(1,e),SM(2,e),SM(3,e),'r*');
    plot3(C(1,:),C(2,:),C(3,:),'b.');    
    plot3(C(1,fidx),C(2,fidx),C(3,fidx),'ro')
    hold off
    view([0 90])
    drawnow
    [eS eC eU(:,e) eE(:,:,e) eL eERR eLAM] = PCA_FIT_FULL_T(C(:,fidx),1);    
end
for e = 2:numel(Xi)
    for d = 1:size(eE,2)
        if eE(:,d,e-1)'*eE(:,d,e) < 0
            eE(:,d,e) = -eE(:,d,e);
        end
    end
end
%% watch paches
close all
for e = 1:size(SM,2)
    M = PCA_BKPROJ_T(SM(:,e),E,U);
    M = reshape(M,[255 3]);
    M = bsxfun(@times,M,sum(M,1).^-1);
    K = (1:255)*M;
    plot(M);
    for k = 1:3
        P(:,:,k) = K(k)*ones(200,200);
    end
    %imshow(P/255,[])
    drawnow
end
%% watch hitograms
close all
for e = 1:size(H,3)
    plot(H(:,:,e))
    drawnow
    pause(.3)
end
%% display the pca values curvilinear
CL = {'m' 'c'};
for e = 1:numel(Xi)
    for d = 1:size(eE,2)
        quiver3(SM(1,e),SM(2,e),SM(3,e),eE(1,d,e),eE(2,d,e),eE(3,d,e),.1,'Color',CL{d})
        drawnow
    end
end

%% view the results of parameterization
% note the first image is as a different resolution
close all
tcH = cell(1,100);
disp = 0;
CR = {};
% make n hood to sample the stem-end which is "always" up
parfor e = 1:100%numel(FileList)
    imgInfo = imfinfo(FileList{e});
    if imgInfo.Width < 4000
        e
        % read the image
        I = imread(FileList{e});
        % get the chip mask
        B = getChipMask_ver0(FileList{e},5000,1000);
        fprintf(['Done with mask \n']);
        % make inner ship mask and then subtract from chip mask to reveal the
        % outter rim
        h = strel('disk',51,0);
        iB = imerode(B,h);        
        fprintf(['Done with erode ring \n']);
        % get the boundaries for the mask
        idB = bwboundaries(iB);
        % create the ring mask
        rM = B - iB;
        % sample the ring mask
        [HI CM BB] = sampleRing(rM,I);
        fprintf(['Done with sample ring \n']);
        % make the total count
        %tcH = cat(3,tcH,HI);
        tcH{e} = HI;

        
        % crop out the images
        for i = 1:size(BB,1)
            CR{e}{i} = imcrop(I,BB(i,:));
        end

        % show the user the chips in order and ask a question



        % get the boundaries for the mask
        dB = bwboundaries(B);
        if disp
            % show the boundaries and the image
            close all
            imshow(I,[]);
            hold on    
            for i = 1:numel(dB)
                plot(dB{i}(:,2),dB{i}(:,1),'r')
            end
            for i = 1:numel(idB)
                plot(idB{i}(:,2),idB{i}(:,1),'g')
            end
            hold off
            title(num2str(e));
           
            
            fidx = find(NM==e);
            for i = 1:size(CM,2)
                text(CM(1,i),CM(2,i),[num2str(kidx(fidx(i))) '--' num2str(loc(fidx(i)))]);
            end
         	drawnow
            
            
        end
        
        
    end
    e
   
end
%% stack and sort the images
STACK = {};
for e = 1:numel(CR)
    for i = 1:numel(CR{e})
        STACK{end+1} = CR{e}{i};
    end
end
STACK = STACK(sidx);
for e = 1:numel(STACK)
    mFileList{e} = ['/home/nate/junk/' num2str(e) '.tif'];
    tmp = imresize(STACK{e},[200 200]);
    imwrite(tmp,mFileList{e});
end
%%
montage(mFileList(1:800),'Size',[20 10],'Indices',1:4:800)



