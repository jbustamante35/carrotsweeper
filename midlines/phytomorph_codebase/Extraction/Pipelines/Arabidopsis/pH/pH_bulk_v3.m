%% 1) look for images
FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/pH-measurements/';
FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/test4/';
FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/';
FilePath = '/home/nate/Downloads/Fusicoccin/';
%FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Gravitropism/Gaby/';
%FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Auxin treatment/Gaby/';
FileList = {};
FileExt = {'tif','TIF'};
verbose = 1;
FileList = sdig(FilePath,FileList,FileExt,verbose);
%% 2) watch movies
SKIP = 10;
for s = 1:numel(FileList)
    for e = 1:SKIP:numel(FileList{s})
        I = imread(FileList{s}{e});
        I = I(:,:,1:3);
        imshow(I,[]);
        drawnow
    end
end
%% 2.5) remove 1122 data
moniker = '1122';
rmidx = [];
for e = 1:numel(FileList)
    [p,n,ext] = fileparts(FileList{e}{1});
    if ~isempty(strfind(p,moniker))
       rmidx = [rmidx e];
    end
end
FileList(rmidx) = [];
%% 2.5) remove non 2015 data
moniker = '2015';
rmidx = [];
for e = 1:numel(FileList)
    [p,n,ext] = fileparts(FileList{e}{1});
    if isempty(strfind(p,moniker))
       rmidx = [rmidx e];
    end
end
FileList(rmidx) = [];
%% 2.5) remove non gravitropism data
moniker = 'Auxin';
moniker = '20150213';
rmidx = [];
for e = 1:numel(FileList)
    [p,n,ext] = fileparts(FileList{e}{1});
    if isempty(strfind(p,moniker))
       rmidx = [rmidx e];
    end
end
FileList(rmidx) = [];
%% find cngc
moniker = 'cngc';
subSet = [];
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e}{1},moniker))
       subSet = [subSet e];
    end
end
%% 2.6) calibrate
close all
phD = readtext('/home/nate/Downloads/pH calibration - 2014-12 (center).csv');
phD = readtext('/home/nate/Downloads/Results_pH calibration for_20150130_pH_horizontal center of FOV.csv');
phD = readtext('/home/nate/Downloads/20150213 pH calibration_vertical stage (central FOV bar).csv');
PH = cell2mat(phD(:,1));
uPH = mean(PH);
RATIO = cell2mat(phD(:,2));
uR = mean(RATIO);
PH = PH - uPH;
RATIO = RATIO - uR;
sigfunc = @(A, x)A(5) + (A(1) ./ (A(2) + exp(-A(3)*(x-A(4)))));
A0 = [1 1 1 1 0]; %// Initial values fed into the iterative algorithm
fit = nlinfit(PH, RATIO, sigfunc, A0);
iPH = linspace(4-uPH,8-uPH,10000);
Y = sigfunc(fit,iPH);
plot(PH+uPH,RATIO+uR,'b*');
hold on
plot(iPH+uPH,Y+uR,'r');
iPH = iPH + uPH;
Y = Y + uR;
%% make image calibration over space
%FilePath = '/home/nate/Downloads/ph_auxin/';
FilePath = '/home/nate/Downloads/ph_gravi/';
FilePath = '/home/nate/Downloads/frameAverageB/';
%FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Gravitropism/Gaby/';
%FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Auxin treatment/Gaby/';
FileList = {};
FileExt = {'tif','TIF'};
verbose = 1;
caliList = gdig(FilePath,FileList,FileExt,verbose);
phVal = [];
for e = 1:numel(caliList)
    [pth nm ext] = fileparts(caliList{e});
    fidx1 = strfind(nm,'pH');
    fidx2 = strfind(nm,' ');
    fidx2(fidx2 < fidx1(1)) = [];
    phVal(e) = str2num(nm(fidx1(1)+3:fidx2(2)-1));
end
[phVal,sidx]= sort(phVal);
caliList = caliList(sidx);
C = [];
scale = 1/8;
for e = 1:numel(caliList)
    tmp = imread(caliList{e});
    tmp = imfilter(tmp,fspecial('average',7),'replicate');
    tmp = imresize(tmp,scale);
    C = cat(4,C,tmp);
    imshow(C(:,:,:,e),[]);
    pause(.1);
    drawnow
end
C = double(C);
uphVal = mean(phVal);
phVal = phVal - uphVal;
fit = [];
% perform fitting
for i = 1:size(C,1)
    parfor j = 1:size(C,2)
        tic
        R = squeeze(C(i,j,:,:));
        R = squeeze(R(2,:).*R(1,:).^-1);
        uR(i,j) = mean(R);
        R = R - uR(i,j);
        sigfunc = @(A, x)A(5) + (A(1) ./ (A(2) + exp(-A(3)*(x-A(4)))));
        A0 = [1 1 1 1 0]; %// Initial values fed into the iterative algorithm
        fit(i,:,j) = nlinfit(phVal, R, sigfunc, A0);
        toc
    end
    i
end
fit = permute(fit,[1 3 2]);
%% view synethic images
close all
ph = phVal + uphVal;
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
cnt = 1;
for p = ph
    Io = double(imread(caliList{cnt}));
    Io = Io(:,:,2).*Io(:,:,1).^-1;
    Io = imfilter(Io,fspecial('average',7),'replicate');
    figure(h3);
    %imshow(Io,[1 3]);
    mesh(Io);
    colorbar
    
    view([0 90]);
    figure(h4);
    hist(Io(:),100);
    R = [];
    for i = 1:size(C,1)
        for j = 1:size(C,2)
            R(i,j) = sigfunc(fit(i,j,:),p-uphVal) + uR(i,j);
        end
    end
    R = imfilter(R,fspecial('average',3),'replicate');
    R = imresize(R,8);
    figure(h3);
    caxis([min(R(:)) max(R(:))])
    figure(h1);
    mesh(R);
    colorbar
    caxis([min(R(:)) max(R(:))])
    view([0 90]);
    %imshow(R,[1 3]);
    title(p);
    drawnow
    figure(h2);
    hist(R(:),100)
    
    waitforbuttonpress
    cnt = cnt + 1;
end
%%  look at single image for ratio on vertical stage
close all
V1 = double(imread('~/Downloads/pH 6 vertcial_0000.tif'));
V = double(imread('~/Downloads/pH 6 vertcial_0007.tif'));
V = double(imread('/home/nate/Downloads/ph_gravi/20150202_pH 5 buffer_verticalstage_t1.tif'));
V1 = double(imread('~/Downloads/20150202_pH 5 buffer_verticalstage_average.tif'));
V1 = double(imread('~/Downloads/20150202_pH 5.2 buffer_verticalstage_16frameaverage.tif'));
V = imfilter(V,fspecial('average',5),'replicate');
V = V(:,:,2).*V(:,:,1).^-1;
V1 = imfilter(V1,fspecial('average',5),'replicate');
V1 = V1(:,:,2).*V1(:,:,1).^-1;
figure;
mesh(V)
view([0 90]);
caxis([1.8 2.2])
figure;
mesh(V1)
view([0 90]);
caxis([2.2 2.6])
%% 3) extract root
disp = 0;
for s = 1:numel(FileList)
    s
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear samp    
    norScale = linspace(0,10,10);
    sampleAlong = 400;
    tdb = {};
    tC = {};
    tN = {};
    tS = {};
    tM = {};
    tT = {};
    tB = {};
    parfor e = 1:numel(FileList{s})
        tic
        fileName = FileList{s}{e};
        tdb{e} = extractBoundary(fileName,11);        
        %tdb{e} = smoothCurve(tdb{e});
        [tC{e}] = getTip(tdb{e},31);
        [tN{e}] = generateNormalField(tdb{e},7);
        [tS{e} tM{e} tT{e} tB{e}] = sampleBoundary(fileName,tdb{e},tN{e},norScale,tC{e},sampleAlong,0);
        toc*numel(FileList{s})/60/12
        s
    end
    dB{s} = tdb;
    tipC{s} = tC;
    NOR{s} = tN;
    samp{s} = tS;
    mask{s} = tM;
    topIdx{s} = tT;
    botIdx{s} = tB;
end
%%
%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display heat map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = [];
    for e = 1:numel(FileList{s})
        S = cat(3,S,samp{s}{e});
    end
    S = permute(S,[1 3 2]);    
    ratio = S(:,:,2).*S(:,:,1).^-1;
    phV = [];
    for e = 1:numel(ratio)
        delta = abs(ratio(e) - Y);
        [~,midx] = min(delta);
        phV(e) = iPH(midx);
    end
    phV = reshape(phV,size(ratio));
    
    RAW{s} = S;
    PHDATA{s} = phV;
    
    if disp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display movie with contour
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nSKIP = 1;
        for e = 1:numel(FileList{s})
            I = imread(FileList{s}{e});
            I = I(:,:,1:3);
            I = permute(I,[2 1 3]);
            imshow(I,[]);
            hold on        
            ED = [];
            for p = 1:size(dB{s}{e}{1},1)
                str = dB{s}{e}{1}(p,:);
                ed = str + norScale(end)*NOR{s}{e}(p,:);
                if mod(p,nSKIP) == 0
                    plot([str(2) ed(2)],[str(1) ed(1)],'r')
                end
                ED = [ED;ed];
            end



            plot(ED(:,2),ED(:,1),'r');
            plot(dB{s}{e}{1}(:,2),dB{s}{e}{1}(:,1),'r');

            plot(ED(topIdx{s}{e}(1):tipC{s}{e},2),ED(topIdx{s}{e}(1):tipC{s}{e},1),'c');
            plot(dB{s}{e}{1}(topIdx{s}{e}(1):tipC{s}{e},2),dB{s}{e}{1}(topIdx{s}{e}(1):tipC{s}{e},1),'c');


            plot(ED(tipC{s}{e}:botIdx{s}{e}(end),2),ED(tipC{s}{e}:botIdx{s}{e}(end),1),'b');
            plot(dB{s}{e}{1}(tipC{s}{e}:botIdx{s}{e}(end),2),dB{s}{e}{1}(tipC{s}{e}:botIdx{s}{e}(end),1),'b');



            plot(dB{s}{e}{1}(tipC{s}{e},2),dB{s}{e}{1}(tipC{s}{e},1),'g*');

            drawnow
            hold off
        end
    end
        
end
%}
%% watch movies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display movie with contour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
nSKIP = 10;
s = 5;
for e = 1:numel(FileList{s})
    I = imread(FileList{s}{e});
    imshow(I,[]);
    hold on        
    ED = [];
    for p = 1:size(dB{s}{e}{1},1)
        str = dB{s}{e}{1}(p,:);
        ed = str + norScale(end)*NOR{s}{e}(p,:);
        if mod(p,nSKIP) == 0
            plot([str(2) ed(2)],[str(1) ed(1)],'r')
        end
        ED = [ED;ed];
    end



    plot(ED(:,2),ED(:,1),'r');
    plot(dB{s}{e}{1}(:,2),dB{s}{e}{1}(:,1),'r');

    plot(ED(topIdx{s}{e}(1):tipC{s}{e},2),ED(topIdx{s}{e}(1):tipC{s}{e},1),'c');
    plot(dB{s}{e}{1}(topIdx{s}{e}(1):tipC{s}{e},2),dB{s}{e}{1}(topIdx{s}{e}(1):tipC{s}{e},1),'c');


    plot(ED(tipC{s}{e}:botIdx{s}{e}(end),2),ED(tipC{s}{e}:botIdx{s}{e}(end),1),'b');
    plot(dB{s}{e}{1}(tipC{s}{e}:botIdx{s}{e}(end),2),dB{s}{e}{1}(tipC{s}{e}:botIdx{s}{e}(end),1),'b');



    plot(dB{s}{e}{1}(tipC{s}{e},2),dB{s}{e}{1}(tipC{s}{e},1),'g*');

    drawnow
    hold off
end
%% display all data, cat WT and cngn, title graphs
close all
WTD = [];
MTD = [];
oPath = '/home/nate/Downloads/phOut/';
mkdir(oPath);
N = 420;
for e = 1:numel(PHDATA)
    figure;
    tmp = PHDATA{e};
    tmp = imfilter(tmp,fspecial('average',11),'replicate');
    mesh(tmp);
    axis([0 size(PHDATA{e},2) 0 800])
    view([0 90]);
    colorbar
    caxis([5 6.6])
    [pth nm ext] = fileparts(FileList{e}{1});
    filename = [oPath nm '.csv'];
    csvwrite(filename,PHDATA{e}(:,1:N));
    saveas(gca,strrep(filename,'csv','tif'));
    %if ~isempty(strfind(FileList{e}{1},'WT')) & ~isempty(strfind(FileList{e}{1},'Auxin'))
    if ~isempty(strfind(FileList{e}{1},'WT'))
        title('WT')
        WTD = cat(3,WTD,PHDATA{e}(:,1:N));
        
    end
    
    %if ~isempty(strfind(FileList{e}{1},'cngc')) & ~isempty(strfind(FileList{e}{1},'Auxin'))
    if ~isempty(strfind(FileList{e}{1},'cngc'))
        title('cngc')
        MTD = cat(3,MTD,PHDATA{e}(:,1:N));
    end
end
%% 
%% delete first data point
WTD(:,:,1) = [];
%% trim first data set
PHDATA{1}(:,1:200) = [];
%%
figure;mesh(mean(WTD,3))
axis([0 N 0 800])
view([0 90]);
colorbar
caxis([5 6.6])
figure;mesh(mean(MTD,3))
axis([0 N 0 800])
view([0 90]);
colorbar
caxis([5 6.6])
%%

%%
uWT = mean(WTD(:,1:80,:),3);
uWT = imfilter(uWT,fspecial('average',5),'replicate');
figure;
mesh(uWT);
axis([0 size(uWT,2) 0 800])
view([0 90]);
colorbar
caxis([4 7]);
uMT = mean(MTD(:,1:80,:),3);
uMT = imfilter(uMT,fspecial('average',5),'replicate');
figure;
mesh(uMT);
axis([0 size(uMT,2) 0 800])
view([0 90]);
colorbar
caxis([4 7]);
figure;
mesh(uWT);
hold on
mesh(uMT);
figure;
mesh(uWT - uMT);