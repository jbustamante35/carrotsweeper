FilePath = '/home/nate/Downloads/CAM/Camera402/';
FileList = {};
FileExt = {'jpg'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
FilePath = '/mnt/snapper/nate/fieldData/';
FileList = {};
FileExt = {'JPG'};
SET = sdig(FilePath,FileList,FileExt,1);
FileList = SET{1};
%%
D = imread(FileList{1});
D = zeros([size(D) numel(FileList)]);
for e = 1:numel(FileList)
    D(:,:,:,e) = imread(FileList{e});
    e
end
%%
for e = 1:numel(FileList)
    imshow(D(:,:,:,e)/255,[]);
    drawnow
end
%% to hsv
for e = 1:size(D,4)
    s(:,:,e)= rgb2hsv_fast(D(:,:,:,e)/255,'','H');
    e
end
%%
for e = 1:size(s,3)
    imshow(s(:,:,e),[]);
    drawnow
end
%% hard filter
for e = 1:size(s,3)
    imshow(s(:,:,e) > .15 & s(:,:,e) < .4,[])
    drawnow
end
%%
G = [];
for e = 1:size(s,3)
    tmp = imresize(s,.5);
    G = [G;tmp(:)];
end
%%
obj = gmdistribution.fit(G(1:100:end),3);
%%
close all
vec = [];
parfor e = 1:size(s,3) 
    %tmp = s(:,:,e);
    tmp = imread(FileList{e});
    tmp = rgb2hsv_fast(tmp,'','H');
    idx = cluster(obj,tmp(:));
    idx = reshape(idx,[size(s,1) size(s,2)]);
    M = idx == 1;
    M = bwareaopen(M,3000);
    M = imclearborder(M);
    vec(:,e) = sum(M,2);
    %{
    imshow(M,[]);
    hold on
    plot(vec(:,e)*.5,1:size(vec,1),'r');
    hold off
    drawnow
    %}
end
%%
vec = [];
parfor e = 1:10%numel(FileList)
    [vec(:,e) tot(e)] = frameMetric_ver0(FileList{e},obj,1);
    
end
%%
dataPath = ['/iplant/home/cmcninch/stationary_images/CAM298%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
FileList = {};
 for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
[FileList] = issueBulkTicket(FileList);
    %% remove if contians thumb.jpg
    rm = [];
    for e = 1:numel(r)
        if ~isempty(strfind(r(e).DATA_NAME,'.thumb.jpg'))
            rm = [rm e];
        end
    end
    r(rm) = [];
%%
func = cFlow('frameMetric_ver0');
func.setMCRversion('v840');
func.setMemory(2000);
numJobs = numel(FileList);
for e = 1:numJobs
    fprintf(['start generating job:' num2str(e) ':' num2str(numJobs) '\n']);
    [m1r{e} m2r{e}] = func(FileList{e},obj,1);
    fprintf(['end generating job:' num2str(e) ':' num2str(numJobs) '\n']);
end
func.submitDag(auth,500,500);
%%
for e = 1:numJobs
    lm1{e} = cFlowLoader(m1r{e});
end