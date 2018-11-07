%% run more data local - scan first
FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/ssu1/Bd21_11_26_2015/';
FileList = {};
FileExt = {'tiff','TIF'};
verbose = 1;
masterSet = sdig(FilePath,FileList,FileExt,verbose);
masterSet{1}(600:end) = [];
%% define output folder
oPath = '/mnt/scratch1/phytoM/outputSSu1_2/';
mkdir(oPath);
%% first part of code snip is for quick local run for grant
% 
% run local
% input list for isolateRoots_overStack
% (stack,outPath,outName,toWrite,NP,NPK,SNIP,disp,numToProcess)
% stack - the file names of the images    
% outPath - the location to write the results to
% outName - the name of the files to save the results as
% toWrite - flag for writing csv files    
% NP - number of points for tip angle measurement
% NPK - number of points to measure kurvature over
% SNIP - number of points to SNIP from midline
% disp - to display
% numImagesToProcess - if less than zero then process all


%%%% run analysis on each set
for s = 1:numel(masterSet)
    tm = clock;
    out{s} = isolateRoots_overStack(masterSet{s},oPath,[num2str(s) '--data'],1,20,20,20,0,-1);
    etm = etime(clock,tm);
end

%% next step is single analysis of grant data
close all
FilePath = '/mnt/scratch1/phytoM/outputSSu1/';
%FilePath = oPath;
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
angleList = {};
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'angle'))
        angleList{end+1} = FileList{e};
    end
end
%angleList([2 6]) = [];
clear adata
for e = 1:numel(angleList)
    adata(e).data = csvread(angleList{e})*180/pi;
    adata(e).name = angleList{e};
end



for e = 1:numel(adata)
    figure;
    X = 1:size(adata(e).data,1);
    X = X*2/60;
    plot(X,adata(e).data)
    
    
    csvwrite('/mnt/spaldingdata/nate/oldData_angle.csv',adata(e).data)
    
end
%close all

%%
% glue
sumData(1).data = [adata(2).data;adata(1).data];
sumData(2).data = [adata(4).data;adata(3).data];

for e = 1:numel(sumData)
    X = 1:size(sumData(e).data,1);
    X = X*2/60;
    figure;
    plot(X,sumData(e).data);  
    
    
end

%csvwrite('/mnt/spaldingdata/nate/d1.csv',sumData(1).data)
%csvwrite('/mnt/spaldingdata/nate/d2.csv',sumData(2).data)
%csvwrite('/mnt/spaldingdata/nate/g1.csv',sumLData(1).data)
%csvwrite('/mnt/spaldingdata/nate/g2.csv',sumLData(2).data)

close all
% fft
SEL = [1 2];
for e = 1:numel(SEL)
    toOp = sumData(SEL(e)).data;
    tmp = imfilter(toOp,fspecial('average',[181 1]),'replicate')
    figure
    plot(tmp)
    hold on
    plot(toOp);
    hold off
    
    f = toOp - tmp;
    figure;plot(f);
    f = bsxfun(@minus,f,mean(f,1));
    ft = abs(fft(f,[],1))
    figure;
    plot(mean(ft,2))
    %waitforbuttonpress
end

growthList = {};
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'growthRate'))
        growthList{end+1} = FileList{e};
    end
end
growthList([2 6]) = [];




clear ldata
for e = 1:numel(growthList)
    ldata(e).data = csvread(growthList{e});
    ldata(e).name = growthList{e};
end


% glue
sumLData(1).data = [ldata(2).data;ldata(1).data];
sumLData(2).data = [ldata(4).data;ldata(3).data];

for e = 1:numel(sumData)
    X = 1:size(sumLData(e).data,1);
    X = X*2/60;
    figure;
    plot(X,sumLData(e).data);    
    figure;
    plot(X,mean(sumLData(e).data,2));
    mean(mean(sumLData(e).data))
end
%%
close all
FilePath = '/mnt/scratch1/phytoM/outputSSu1/';
%FilePath = oPath;
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);

matList = {};
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'mat'))
        matList{end+1} = FileList{e};
    end
end
%matList([2 6]) = [];
%%
A = load(matList{5});

para = {[5]};
NP = 200;
SNIP = 10;
K = [];
for m = 1:(numel(A.out{1}.midlines)-1)
    for t = 1:numel(A.out)
        curve = A.out{t}.midlines(m).data;
        o = cwtK_filter(curve',para);
        K(:,t,m) = o.K(SNIP:NP);
    end
end
%%
close all
mesh(K(:,:,1))
%%
figure
SEL = 2;


for t = 1:numel(A.out)
    plot(A.out{t}.midlines(SEL).data(1,1:10),-A.out{t}.midlines(SEL).data(2,1:10),'r')
    hold on
    plot(A.out{t}.contours(SEL).data(1,:),-A.out{t}.contours(SEL).data(2,:),'b')
    axis equal
    drawnow
    hold off
end




