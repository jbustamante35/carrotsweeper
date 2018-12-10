%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scan for new images - STEP 1
inFilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
%inFilePath = '/mnt/piranhapuffer2/Bessie/MaizeGravitropism/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(inFilePath,FileList,FileExt,verbose);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run local - STEP 2.1
basePath = '/mnt/spaldingdata/nate/mJUNK/';
mkdir(basePath);
for e = 1:numel(SET)
    % generate save name
    [pth nm ext] = fileparts(SET{e}{1});
    
    
    pth = strrep(pth,filesep,'SLASH');
    saveName = [pth '.mat'];
    saveName = strrep(saveName,' ','SPACE');
    
    
    tic
    cropMaizeImage(SET{e},saveName,basePath,0,1,'/mnt/spaldingdata/nate/labelData3.mat');   
    toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run condor - STEP 2.2
toCompile = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if toCompile
    compile_directory = '/mnt/scratch1/phytoM/flashProjects/maize/';
    CMD = ['mcc -d ' compile_directory ' -m -v -R -singleCompThread cropMaizeImage.m -a /mnt/scratch1/phytoM/phytoSS/matlab/libs/matlab_libs/myHMM/ '];
    eval(CMD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% condor launch - START
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look for image data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%% scan for new images
inFilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
%inFilePath = '/mnt/piranhapuffer2/Bessie/MaizeGravitropism/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(inFilePath,FileList,FileExt,verbose);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove non 61
ridx = [];
for e = 1:numel(SET)
    if numel(SET{e}) < 61
        ridx(e) = 1;
    else
        ridx(e) = 0;
    end
end
SET(find(ridx)) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the dag 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate dag
tmpFileLocation = '/mnt/scratch1/phytoM/flashProjects/maize/';
dag = epfod();
dag.setFunctionName('cropMaizeImage');
%dateString = strrep(strrep(strrep(datestr(clock, 'yyyy-mm-dd HH:MM:SS'),' ','_'),':','_'),'-','_');
outputLocation = ['/mnt/snapper/share/nate/myDev/maizeWhole_mini8/'];
%outputLocation = ['/mnt/snapper/nate/myDev/maizeWhole_mini5/'];
mkdir(outputLocation);
dag.setOutputLocation(outputLocation);
dag.setTempFilesLocation(tmpFileLocation);

numJobs = numel(SET);
toRun = {};
%numJobs = 1;
for e = 1:numJobs
    
    %%%%%%%%%%%%%%%%%%%%%
    % file(s) into parts
    clear condorLocalNames
    for f = 1:numel(SET{e})        
        [pth,condorLocalNames{f},ex] = fileparts(SET{e}{f});
        condorLocalNames{f} = [condorLocalNames{f} ex];
    end
    condorLocalNames = cellStringConvert(condorLocalNames);
    
    %%%%%%%%%%%%%%%%%%%%%
    % generate output name
    pth = strrep(pth,filesep,'SLASH');
    saveName = [pth '.mat'];
    saveName = strrep(saveName,' ','SPACE');

    if ~exist([strrep(outputLocation,'share/','') 'output/' saveName])
        toRun{end+1} = SET{e};
        SET{e}{1}
    end
    %{
    %%%%%%%%%%%%%%%%%%%%%
    % create job
    job = cJob();
    job.setOSG(1);
    job.setFlock(1);
    job.setTempFilesLocation(tmpFileLocation);
    job.setFunctionName('cropMaizeImage');    
    job.setNumberofArgs(6);
    
    %%%%%%%%%%%%%%%%%%%%%
    % add image file(s) to job
    for f = 1:numel(SET{e})
        job.addFile(SET{e}{f});
    end
    % add hmm and pcaData file to the job
    job.addFile('/mnt/spaldingdata/nate/labelData3.mat');
    
    
    %%%%%%%%%%%%%%%%%%%%%
    % set arguments
    job.setArgument(condorLocalNames,1);
    job.setArgument(saveName,2);
    job.setArgument('./output/',3);
    job.setArgument('0',4);
    job.setArgument('1',5);
    job.setArgument('labelData3.mat',6);
    
    % add job to dag
    dag.addJob(job);
    job.generate_submitFilesForDag();
    %}
end



dag.renderDagFile();
scpList = dag.generate_scpFileList();
dirCMD_logs_out = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stdout/'''];
dirCMD_logs_err = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stderr/'''];
dirCMD_output = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/output/'''];
[status result] = system(strrep(dirCMD_logs_out,'#directory#',dag.jobFunction));
[status result] = system(strrep(dirCMD_logs_err,'#directory#',dag.jobFunction));
[status result] = system(strrep(dirCMD_output,'#directory#',dag.jobFunction));
dirCMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir /home/nate/condorFunctions/#directory#/'''];
[status result] = system(strrep(dirCMD,'#directory#',dag.jobFunction));
CMD = 'scp -P 50118 #srcfile# nate@128.104.98.118:/home/nate/condorFunctions/#directory#/#desfile#';
CMD = strrep(CMD,'#directory#',dag.jobFunction);
for f = 1:numel(scpList)
    [pth nm ext] = fileparts(scpList{f});
    tCMD = strrep(CMD,'#desfile#',[nm ext]);
    tCMD = strrep(tCMD,'#srcfile#',scpList{f});
    [status result] = system(tCMD);
end

% submit the job dag
dagName = dag.generate_dagName();
CMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'cd /home/nate/condorFunctions/#directory#/; condor_submit_dag -maxidle 75 -maxjobs 75 -maxpost 50 ' dagName ''''];
CMD = strrep(CMD,'#directory#',dag.jobFunction);
system(CMD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% condor launch - END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run local - STEP 2.1
basePath = '/mnt/spaldingdata/nate/mJUNK/';
mkdir(basePath);
for e = 1:numel(toRun)
    % generate save name
    [pth nm ext] = fileparts(toRun{e}{1});
    
    
    pth = strrep(pth,filesep,'SLASH');
    saveName = [pth '.mat'];
    saveName = strrep(saveName,' ','SPACE');
    
    
    tic
    cropMaizeImage(toRun{e},saveName,basePath,0,1,'/mnt/spaldingdata/nate/labelData3.mat');   
    toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load from condor bulk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for generating <HMM>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inFilePath = ['/mnt/snapper/nate/myDev/maizeWhole/'];
FileList = {};
FileExt = {'mat'};
verbose = 1;
mFileList = gdig(inFilePath,FileList,FileExt,verbose);
close all
dispi = 1;
frDS = 30;
curveData = [];
dcurveData = [];
patchData = [];
nameData = [];
frameData = [];
cnt = 1;
MAXLOAD = 500;
for m = 1:numel(mFileList)
    load(mFileList{m},'numberFrames','numberSeedlings');
    for s = 1:numberSeedlings
        for tm = 1:frDS:numberFrames
            try
                curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                frameString = ['seedling' num2str(s) 'frame' num2str(tm)];
                f = load(mFileList{m},frameString,curveString);


                sz = size(f.(curveString).S);
                patchData = [patchData;reshape(f.(curveString).S,[sz(1)*sz(2) sz(3)])'];
                curveData = [curveData;f.(curveString).data'];
                nameData = [nameData;cnt*ones(sz(3),1)];
                dcurveData = [dcurveData;bsxfun(@minus,f.(curveString).data,mean(f.(curveString).data,2))'];
                frameData = cat(3,frameData,f.(curveString).E);
                cnt = cnt + 1;


                if dispi
                    imshow(f.(frameString),[]);
                    hold on
                    plot(f.(curveString).data(1,:),f.(curveString).data(2,:),'r');
                    hold off
                    drawnow
                end
            
                fprintf(['Done dataset ' num2str(m) ':' num2str(numel(mFileList)) ':' num2str(cnt) '\n']);
                
            catch ME
            
            
            end
            
            
        end
       
    end
    if cnt >= MAXLOAD
        break
    end
end
%% PCA on patch data
[pS pC pcaData.pU pcaData.pE pL pERR pLAM] = PCA_FIT_FULL(patchData,2);
%% segment on patches via gmm
k = 4;
learnDims = 2;
opt = statset('display','iter','MaxIter',500);
obj = gmdistribution.fit(pC(:,1:learnDims),k,'Options',opt);
%% display segmentation from gmm
close all
UQ = unique(nameData);
h1 = figure;
CL = {'r' 'g' 'b' 'k' 'y' 'm'};
for u = 1:numel(UQ)
    figure(h1);    
    fidx = find(nameData==UQ(u));
    Psig = pC(fidx,1:learnDims);    
    
    
    y = [];
    for g = 1:size(obj.mu,1)
        y(:,g) = mvnpdf(Psig,obj.mu(g,:),obj.Sigma(:,:,g));
    end
    [~, gidx] = max(y,[],2);
    
       
    tmpC = curveData(fidx,:);
    plot(tmpC(:,1),tmpC(:,2));
    hold on
    % get tangent and normal space
    E = frameData(:,:,fidx);    
    
    quiver(tmpC(:,1),tmpC(:,2),squeeze(E(1,2,:)),squeeze(E(2,2,:)));
    for g = 1:size(obj.mu,1)
        fidx = find(gidx==g);
        plot(tmpC(fidx,1),tmpC(fidx,2),[CL{g} '*']);
    end
    axis equal
    hold off
    waitforbuttonpress;
end
%% segment from gmm and ...
y = [];
for g = 1:size(obj.mu,1)
    y(:,g) = mvnpdf(pC(:,1:learnDims),obj.mu(g,:),obj.Sigma(:,:,g));
end
[~, gidx_master_Y] = max(y,[],2);
%% get parameters for distributions
tipGrp = 4;
rootGrps = 3;
transGrps = 1;
kernelGrp = 2;

tip_mean_patch = mean(pC(gidx_master_Y==tipGrp,1:learnDims));
tip_cov_patch = cov(pC(gidx_master_Y==tipGrp,1:learnDims));
root_mean_patch = mean(pC(gidx_master_Y==rootGrps,1:learnDims));
root_cov_patch = cov(pC(gidx_master_Y==rootGrps,1:learnDims));
trans_mean_patch = mean(pC(gidx_master_Y==transGrps,1:learnDims));
trans_cov_patch = cov(pC(gidx_master_Y==transGrps,1:learnDims));
kernel_mean_patch = mean(pC(gidx_master_Y==kernelGrp,1:learnDims));
kernel_cov_patch = cov(pC(gidx_master_Y==kernelGrp,1:learnDims));
close all
fidx = find(gidx_master_Y==tipGrp);
kidx = kmeans(dcurveData(fidx,:),2);
plot(dcurveData(fidx(kidx==1),1),dcurveData(fidx(kidx==1),2),'b.');
hold on
plot(dcurveData(fidx(kidx==2),1),dcurveData(fidx(kidx==2),2),'r.');
tip_mean_cord = mean(dcurveData(fidx(kidx==1),:));
tip_cov_cord = cov(dcurveData(gidx_master_Y==tipGrp,:));
root_mean_cord = mean(dcurveData(gidx_master_Y==rootGrps,:));
root_cov_cord = cov(dcurveData(gidx_master_Y==rootGrps,:));
trans_mean_cord = mean(dcurveData(gidx_master_Y==transGrps,:));
trans_cov_cord = cov(dcurveData(gidx_master_Y==transGrps,:));
kernel_mean_cord = mean(dcurveData(gidx_master_Y==kernelGrp,:));
kernel_cov_cord = cov(dcurveData(gidx_master_Y==kernelGrp,:));
%% generate HMM
tip_node = hmm_node('tip');
root_node_upper = hmm_node('root_upper');
t1_node = hmm_node('ut');
kernel_node = hmm_node('kernel_node');
t2_node = hmm_node('lt');
root_node_lower = hmm_node('root_lower');
kernel_node_end = hmm_node('e_kerel_node');



tf1 = heavisideTransitionFunction(20,@(x,y)lt(x,y));
tf2 = heavisideTransitionFunction(20,@(x,y)ge(x,y));


t1 = heavisideTransitionFunction(20,@(x,y)lt(x,y));
t2 = heavisideTransitionFunction(20,@(x,y)ge(x,y));
kernel_to_transZone = constantTransitionFunction(.1);
kernel_to_kernel = constantTransitionFunction(.9);
kernelEnd_to_kernelEnd = constantTransitionFunction(1);
root_to_root = constantTransitionFunction(.75);
root_to_outofToot = constantTransitionFunction(.25);
tip_to_tip = constantTransitionFunction(.75);
tip_to_lower_root = constantTransitionFunction(.25);

kernel_node.attachNode(kernel_node,kernel_to_kernel);
kernel_node.attachNode(t1_node,kernel_to_transZone);

t1_node.attachNode(t1_node,tf1);
t1_node.attachNode(root_node_upper,tf2);

root_node_upper.attachNode(root_node_upper,root_to_root);
root_node_upper.attachNode(tip_node,root_to_outofToot);

%tip_node.attachNode(tip_node,t1);
%tip_node.attachNode(root_node_lower,t2);

tip_node.attachNode(tip_node,tip_to_tip);
tip_node.attachNode(root_node_lower,tip_to_lower_root);

root_node_lower.attachNode(root_node_lower,root_to_root);
root_node_lower.attachNode(t2_node,root_to_outofToot);

t2_node.attachNode(t2_node,tf1);
t2_node.attachNode(kernel_node_end,tf2);

kernel_node_end.attachNode(kernel_node_end,kernelEnd_to_kernelEnd);

% attach distributions for patches
tipD = myProb(tip_mean_patch,tip_cov_patch);
tip_node.attachDistribution(tipD,1);
rootD = myProb(root_mean_patch,root_cov_patch);
root_node_lower.attachDistribution(rootD,1);
root_node_upper.attachDistribution(rootD,1);
transD = myProb(trans_mean_patch,trans_cov_patch);
t1_node.attachDistribution(transD,1);
t2_node.attachDistribution(transD,1);
kernelD = myProb(kernel_mean_patch,kernel_cov_patch);
kernel_node.attachDistribution(kernelD,1);
kernel_node_end.attachDistribution(kernelD,1);
% attach distributions for cordinates
tipD_cord = myProb(tip_mean_cord,tip_cov_cord);
tip_node.attachDistribution(tipD_cord,2);
rootD_cord = myProb(root_mean_cord,root_cov_cord);
root_node_lower.attachDistribution(rootD_cord,2);
root_node_upper.attachDistribution(rootD_cord,2);
transD_cord = myProb(trans_mean_cord,trans_cov_cord);
t1_node.attachDistribution(transD_cord,2);
t2_node.attachDistribution(transD_cord,2);
kernelD_cord = myProb(kernel_mean_cord,kernel_cov_cord);
kernel_node.attachDistribution(kernelD_cord,2);
kernel_node_end.attachDistribution(kernelD_cord,2);
% create hmm
hmm = my_hmm();
hmm.addNode(kernel_node);
hmm.addNode(t1_node);
hmm.addNode(root_node_upper);
hmm.addNode(tip_node);
hmm.addNode(root_node_lower);
hmm.addNode(t2_node);
hmm.addNode(kernel_node_end);
hmm.dn = [1 1];
%% clear class hmm and node_hnm
clear tip_node root_node_upper t1_node kernel_node t2_node root_node_lower kernel_node_end hmm kernel_to_transZone kernel_to_kernel kernelEnd_to_kernelEnd root_to_root root_to_outofToot
clear tf1 tf2
clear class my_hmm 
clear class hmm_node
clear class myProb
clear class myTransitionFunction
clear class constantTransitionFunction
clear class heavisideTransitionFunction
%% simulate
[n1 n2] = ndgrid(-300:300,-300:300);
N = [n1(:) n2(:)];
mPDF = zeros(size(n1));
for s = 1:numel(hmm.NodeList)
    PDF = hmm.NodeList{s}.Distribution{2}.getProb(N',1);
    PDF = reshape(PDF,size(n1));
    PDF = bindVec(PDF);
    imshow(PDF',[]);
    title(hmm.NodeList{s}.StateName);
    waitforbuttonpress
    mPDF = mPDF + PDF;
end

imshow(mPDF',[]);
%% display segmentation from hmm
close all
UQ = unique(nameData);
h1 = figure;
CL = {'r' 'b' 'g' 'k' 'y' 'c' 'm' 'c'};
recalc = 0;

for u = 1:100%numel(UQ)
    figure(h1);    
    fidx = find(nameData==UQ(u));
    Psig_patch = pC(fidx,1:learnDims);    
    Psig_cord = dcurveData(fidx,1:learnDims);
    Psig_tot = [Psig_patch';Psig_cord'];
    observation_labels = [ones(size(Psig_patch,2),1);2*ones(size(Psig_cord,2),1)];
    
    
    if recalc
        if u == 1
            gidx = {};
        end
        gidx{u} = hmm.Viterbi(Psig_tot,observation_labels,1);
    end
       
    tmpC = curveData(fidx,:);
    plot(tmpC(:,1),tmpC(:,2));    
    % get tangent and normal space
    E = frameData(:,:,fidx);    
    hold on
    
    
    quiver(tmpC(:,1),tmpC(:,2),squeeze(E(1,2,:)),squeeze(E(2,2,:)),'b');
    UQg = unique(gidx{u});
    for g = 1:numel(UQg)
        fidx = find(gidx{u}==UQg(g));
        plot(tmpC(fidx,1),tmpC(fidx,2),[CL{UQg(g)} '*']);
    end
    %plot(tmpC(1,1),tmpC(1,2),'ro');
    %plot(tmpC(2,1),tmpC(2,2),'go');
    axis equal
    hold off
    drawnow
    pause(.1)
    %waitforbuttonpress;
end
%% loop train
for e = 1:2

    e
    
    tic
    close all
    UQ = unique(nameData);
    UD = {};
    gidx = {};
    parfor u = 1:numel(UQ)    
        fidx = find(nameData==UQ(u));
        Psig_patch = pC(fidx,1:learnDims);    
        Psig_cord = dcurveData(fidx,1:learnDims);
        Psig_tot = [Psig_patch';Psig_cord'];
        observation_labels = [ones(size(Psig_patch,2),1);2*ones(size(Psig_cord,2),1)];
        gidx{u} = hmm.Viterbi(Psig_tot,observation_labels,1);
        UD{u} = [gidx{u};Psig_tot];
        u
    end
    toc
    
    
    
    observation_labels = [ones(learnDims,1);2*ones(size(dcurveData,2),1)];
    hmm.update(UD,observation_labels);
    
end
%% save hm and pca data
save('/mnt/spaldingdata/nate/labelData3.mat','hmm','pcaData');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for generating </HMM>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% label curves TEST
for m = 1:numel(mFileList)
    load(mFileList{m},'numberFrames','numberSeedlings');
    for s = 1:numberSeedlings
        for tm = 1:numberFrames
            try
                curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                frameString = ['seedling' num2str(s) 'frame' num2str(tm)];
                f = load(mFileList{m},frameString,curveString);
                [labels{s,tm} tip{s,tm} tipVec{s,tm} tipAngle(s,tm)] = labelRoots(f.(curveString),pcaData,hmm,@(x)min(x),1,2,25);
                fprintf(['Done dataset ' num2str(tm) ':' num2str(numberFrames) ':' num2str(s) ':' num2str(numberSeedlings) '\n']);
                
            catch ME
            
            
            end
            
            
        end
       
    end
    
    for s = 1:numberSeedlings
     for tm = 1:numberFrames
            try
                curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                frameString = ['seedling' num2str(s) 'frame' num2str(tm)];
                f = load(mFileList{m},frameString,curveString);
                imshow(f.(frameString),[]);hold on;



                tmpC = f.(curveString).data';
                plot(tmpC(:,1),tmpC(:,2));    
                % get tangent and normal space
                E = f.(curveString).E;    
                hold on


                quiver(tmpC(:,1),tmpC(:,2),squeeze(E(1,2,:)),squeeze(E(2,2,:)),'b');
                quiver(tmpC(tip{tm},1),tmpC(tip{tm},2),tipVec{tm}(1),tipVec{tm}(2),100,'Color','r');
                UQg = unique(labels{tm});
                for g = 1:numel(UQg)
                    fidx = find(labels{tm}==UQg(g));
                    plot(tmpC(fidx,1),tmpC(fidx,2),[CL{UQg(g)} '*']);
                end


                drawnow

            catch ME


            end

     end
    end

end
%% load from condor mini - post hmm - startLOAD here
close all
inFilePath = ['/mnt/snapper/nate/myDev/maizeWhole_mini4/output/'];
FileList = {};
FileExt = {'mat'};
verbose = 1;
mFileList = gdig(inFilePath,FileList,FileExt,verbose);

frDS = 1;
dispi = 0;
fdispi = 0;
CL = {'r' 'b' 'g' 'k' 'y' 'c' 'm' 'c'};
close all
MasterTipAngle = [];
MasterGrowthSur = [];
masterName = {};
setupName = {};
plateName = {};
allLost = {};
seedlingNumber = [];
matFile = {};
iPath = {};
fprintf(['Start\n']);
for m = 1:numel(mFileList)
    m
    numel(mFileList)
    load(mFileList{m},'numberFrames','numberSeedlings','CROPBOX');
    tipAngle = [];
    
    try
        % regenerate file name and extract wells
        %%%%%%%%%%%%%%%%%%%%%
        names = {};
        tmpPlateName = '';
        [pth,nm,ext] = fileparts(mFileList{m});
        % replace mat file name with space
        imagePath = strrep(nm,'SPACE',' ');
        % replace mat file name with slash
        imagePath = strrep(imagePath,'SLASH',filesep); 
        imagePath = strrep(imagePath,[filesep ')'],'-)');
        imagePath = [imagePath filesep];
        fidx = findstr(imagePath,filesep);
        clip = imagePath(fidx(end-1)+1:end-1);
        SETUP = clip(1:2);
        clip(1:3) = [];
        clip = ['_' clip '_'];
        fidx = findstr(clip,'_');
        for e = 2:numel(fidx)-1
            tmpNM = clip(fidx(e)+1:fidx(e+1)-1);
            if numel(tmpNM) == 2
                names{end+1} = tmpNM;
            else
                tmpPlateName = tmpNM;
            end
        end
        tmpPlateName = clip(fidx(1)+1:fidx(2)-1);


        if fdispi 
            % load the first image
            [pth,nm,ext] = fileparts(mFileList{m});
            imagePath = strrep(nm,'SPACE',' ');
            imagePath = strrep(imagePath,'-',filesep); 
            imagePath = strrep(imagePath,[filesep ')'],'-)');
            imagePath = [imagePath filesep];    
            FileList = {};
            FileExt = {'TIF','tif'};
            verbose = 1;
            iFileList = gdig(imagePath,FileList,FileExt,verbose);
            I = imread(iFileList{1});
            imshow(I,[]);
            drawnow
            waitforbuttonpress
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if the number of names is empty or there is not a match for the
        % number of names then all is lost
        if numel(names) ~= numberSeedlings | isempty(names)
            allLost{end+1} = mFileList{m};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load the information about each seedling and over each frame
        if numel(names) == numberSeedlings & ~isempty(names)
            for s = 1:numberSeedlings
                try
                    % init the load strings for the data
                    curveLabelsString = ['seedling' num2str(s) 'curveLabels'];
                    tipAngleString = ['seedling' num2str(s) 'angle'];
                    tipIndexString = ['seedling' num2str(s) 'tipIndex'];
                    tipVecString = ['seedling' num2str(s) 'tipVec'];        
                    % call the load and pout data into fW
                    fW = load(mFileList{m},curveLabelsString,tipAngleString,tipIndexString,tipVecString);
                    % stack the tip angle 
                    tipAngle = [tipAngle;fW.(tipAngleString)];
                    masterName = {masterName{:} names{s}};
                    setupName = {setupName{:} SETUP};
                    plateName = {plateName{:} tmpPlateName};
                    seedlingNumber = [seedlingNumber s];
                    matFile{end+1} = mFileList{m};
                    iPath{end+1} = imagePath;
                    %{
                    % load the curve for obtaining the growthrate
                    tipPath = [];
                    for tm = 1:frDS:numberFrames
                        curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                        tmpData = load(mFileList{m},curveString);
                        tipPath = [tipPath tmpData.(curveString).data(:,fW.(tipIndexString)(tm))];
                    end
                    dP = diff(tipPath,1,2);
                    dP = cumsum([0 sum(dP.*dP,1).^.5]);
                    %}
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % show the movie if dispi is on
                    if dispi
                        for tm = 1:frDS:numberFrames
                            try
                                curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                                f = load(mFileList{m},curveString);
                                f.(curveString).data = bsxfun(@plus,[200 50]'+CROPBOX{s}(1:2)',f.(curveString).data);
                                if dispi
                                    I = imread(iFileList{tm});
                                    imshow(I,[]);
                                    hold on
                                    plot(f.(curveString).data(1,:),f.(curveString).data(2,:),'r');
                                    hold on                    
                                    UQg = unique(fW.(curveLabelsString){tm});
                                    for g = 1:numel(UQg)
                                        fidx = find(fW.(curveLabelsString){tm}==UQg(g));
                                        plot(f.(curveString).data(1,fidx),f.(curveString).data(2,fidx),[CL{UQg(g)} '*']);
                                    end
                                    hold off
                                    axis equal
                                    drawnow
                                end
                                fprintf(['Done dataset ' num2str(m) ':' num2str(numel(mFileList)) ':' num2str(cnt) '\n']);
                            catch ME

                            end

                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % end dispi
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % catch error and display loading error
                catch ME
                    
                    mFileList{s};
                    disp(ME.getReport);
                end
                MasterGrowthSur = [MasterGrowthSur;dP];
            end
            
            plot(180/pi*tipAngle');
            axis([0 61 -30 90]);
            title(num2str(numberSeedlings));
            drawnow
            
            %waitforbuttonpress;
            MasterTipAngle = [MasterTipAngle;tipAngle];
           
        end
    catch ME
        ME
    end
end
%% clean tip angle data
close all
plot(MasterTipAngle');
% filter
dA = abs(diff(MasterTipAngle,1,2));
rmidx = find(any(dA > 15*pi/180,2) | (MasterTipAngle(:,1) < -45*pi/180') | all(MasterTipAngle==0,2) | (MasterTipAngle(:,1) > 60*pi/180'));
%rmidx = find(any(dA > 10*pi/180,2) | all(MasterTipAngle==0,2));
per = (size(MasterTipAngle,1)-numel(rmidx))/size(MasterTipAngle,1);

fMA = MasterTipAngle;
fMG = MasterGrowthSur;


rmData = MasterTipAngle(rmidx,:);
badSeedlingNumber = seedlingNumber(rmidx);
badMat = matFile(rmidx);


fMA(rmidx,:) = [];
%fMG(rmidx,:) = [];
masterName(rmidx) = [];
setupName(rmidx) = [];
plateName(rmidx) = [];
matFile(rmidx) = [];
iPath(rmidx) = [];
size(plateName)
size(setupName)
size(masterName)
size(fMA)
%{
figure;
for e = 1:size(rmData,1)
    plot(rmData(e,:)*180/pi);
    axis([0 61 -45 90])
    drawnow
end
%}
figure;
plot(fMA'*180/pi);
title([num2str(per) '-' num2str(size(fMA,1)) '--' num2str(size(MasterTipAngle,1))]);
%% spool out for Takeshi
oPath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/fullExtract_14.10.14_growthRate/';
mkdir(oPath);
UQ = unique(matFile);
for u = 1:numel(UQ)
    fidx = find(strcmp(matFile,UQ{u}));
    [pth,nm,ext] = fileparts(UQ{u});
    fileName = [oPath nm '.csv'];
    fileName = strrep(fileName,'SLASH','_');
    fileName = strrep(fileName,'SPACE',' ');
    %tmpTipdata = fMA(fidx,:);
    tmpTipdata = fMG(fidx,:);
    csvwrite(fileName,tmpTipdata');
    %plot(tmpTipdata');
    %drawnow
end
%% pour data into csv
DATA_BIG = {};
cnt = 1;
for e = 1:numel(plateName) 
    fidx = strfind(iPath{e},'Moon');    
    if ~isempty(fidx)
        DATA_BIG{cnt,1} = iPath{e}(42:end);
        for q = 1:61
            DATA_BIG{cnt,q+1} = MasterTipAngle(e,q);
        end
        cnt = cnt + 1;
    end
end
cell2csv('/mnt/spaldingdata/nate/vis.csv',DATA_BIG)
%% match spec to gravi
load('/mnt/scratch1/phytoM/flashProjects/maize/specKey.mat','k','g','p');
%plateName = plateNameBK;
load('/mnt/scratch1/phytoM/flashProjects/maize/specKey2.mat','kernelVec','genoVec','popVec','kernelID');
k = kernelVec;
g = genoVec;
p = popVec;
keySet = k.keySet;
itr = keySet.iterator();
tipData = [];
grData = [];
specData = [];
oldTip = [];
prediction = [];
MIS = {};
MISW = {};
MISK = {};
genoType = {};
pl = {};
wn = {};
POP = {};
kerID = [];
kT = 1;
while itr.hasNext()
    kT = kT + 1;
    curKey = itr.next();
    fidx = strfind(curKey,'*');
    curPlateName = curKey(1:fidx-1);
    curPlateName = strrep(curPlateName,'WISN10-','WISN10_');
    curPlateName = strrep(curPlateName,'WISN11-','WISN11_');
    curWellName = curKey(fidx+1:end);
    curWellName = [curWellName(2) curWellName(4)];
    fidx = find(strcmp(plateName,curPlateName) & strcmp(masterName,curWellName));
    if numel(fidx) == 0
        MIS{end+1} = curPlateName;
        MISW{end+1} = curWellName;
        MISK{end+1} = curKey;
        %break
    end
    if numel(fidx) == 1
        try
            tmpData = k.get(curKey)';
            tipData = [tipData;fMA(fidx,:)];
            grData = [grData;fMG(fidx,:)];
            %oldTip = [oldTip;tmpData(780:end)];        
            %specData = [specData;tmpData(1:end-7)];
            specData = [specData;tmpData];
            prediction = [prediction;tmpData(end-6:end)];
            genoType{end+1} = g.get(curKey);
            pl{end+1} = curPlateName;
            wn{end+1} = curWellName;
            POP{end+1} = p.get(curKey);
            kerID = [kerID;kernelID.get([curKey '--' p.get(curKey)])];
        catch
        end
    end
    size(specData)
end
%% test hashmap
%SS = csvread(['/mnt/scratch1/phytoM/flashProjects/maize/shape_spec.csv']);
import java.util.HashMap;
testVec = HashMap();
for e = 1:size(SS,1)
    testVec.put('in',SS(e,:));
    PL = testVec.get('in')';    
    if all(PL == SS(e,:))
        report = 1
    else
        
    end
end
%% match spec to shape
%SS = csvread(['/mnt/scratch1/phytoM/flashProjects/maize/shape_spec.csv']);
%fullD = load(['/mnt/scratch1/phytoM/flashProjects/maize/shape_spec.csv'],'fullD');
%SS = round(SS*1000)/1000;
%toRM = ones(size(specData,1),1);
mySH = [];
subGravi = [];
subLAT = [];
for e = 1:size(specData,1)    
    %toMatch = specData(e,:);    
    toMatch = kerID(e);
    %delta = bsxfun(@eq,SS(:,1:end-4),toMatch);
    %delta = all(delta,2);
    %midx = find(delta);
    midx = find(fullD(:,end)==toMatch);
    
    
    plot(specData(e,:),'b');
    hold on
    plot(SS(midx,4:end-4),'r')
    hold off
    drawnow

    if ~isempty(midx)
        mySH = [mySH ; SS(midx,end-3:end)];
        subGravi = [subGravi ; tipData(e,:)];
        subLAT = [subLAT ; LAT(e,:)];
    end
        %{
        
    %}
        
    e
    size(specData,1);
end
%% look at kmean of init tip angle
close all
clear gS gC gU gE gL gERR gLAM
kidx = kmeans(tipData(:,1),5);
UQ = unique(kidx);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(tipData,3);
nC = [];
% order the groups via mean init value
grU = [];
for u = 1:numel(UQ)
    fidx = kidx == UQ(u);
    subSet = tipData(fidx,1);
    grU(u) = mean(subSet);    
end
[J,nkidxL] = sort(grU);
nkidx = zeros(size(kidx));
for u = 1:numel(UQ)
    fidx = find(kidx == nkidxL(u));
    nkidx(fidx) = u;
end

kidx = nkidx;

for u = 1:numel(UQ)
    fidx = kidx == UQ(u);
    subSet = tipData(fidx,:);
    [gS gC{u} gU(:,u) gE(:,:,u) gL gERR gLAM] = PCA_FIT_FULL(tC(fidx,:),3);
    nC(fidx,:) = gC{u};
    U = mean(tipData(fidx,:),1);
    SE = std(tipData(fidx,:),1,1).*sum(fidx).^-.5;
    errorbar(U,SE)
    hold all
end

gU = gU';
figure;
for u = 1:numel(UQ)
    fidx = kidx == UQ(u);
    plot3(tC(fidx,1),tC(fidx,2),tC(fidx,3),'.','MarkerSize',2);
    hold all
    quiver3(gU(u,1),gU(u,2),gU(u,3),gE(1,1,u),gE(2,1,u),gE(3,1,u),'r');
    quiver3(gU(u,1),gU(u,2),gU(u,3),gE(1,2,u),gE(2,2,u),gE(3,2,u),'g');
    quiver3(gU(u,1),gU(u,2),gU(u,3),gE(1,3,u),gE(2,3,u),gE(3,3,u),'b');
end
dL = diff(gU,1,1);
dL = sum(dL.*dL,2).^.5;
L = cumsum([0;dL]);
sp = spap2(1,5,L',gU');

along = linspace(-.1*max(L),max(L)+.1*max(L),100);
BB = fnval(sp,along)';
BB = interp1(L,gU,along,'pchip','extrap');
plot3(BB(:,1),BB(:,2),BB(:,3),'k','LineWidth',3);
% attach to each point along curve and build up manifold
for e = 1:size(tC,1)
    delta = bsxfun(@minus,BB,tC(e,:));
    delta = sum(delta.*delta,2);
    [~,aidx(e)] = min(delta);
    lav(e) = along(aidx(e));
end
dTH = .6;

dBB = gradient(BB');
dBB2 = gradient(gradient(BB'));
dBBlen2 = sum(dBB2.*dBB2,1).^-.5;
dBBlen = sum(dBB.*dBB,1).^-.5;
dBB = bsxfun(@times,dBB,dBBlen);
dBB2 = bsxfun(@times,dBB2,dBBlen2);
alongTAN = bsxfun(@times,dBB,dot(dBB2,dBB,1));
NOR = dBB2 - alongTAN;
NORlen = sum(NOR.*NOR,1).^-.5;
NOR = bsxfun(@times,NOR,NORlen);
% build out frames
C = cross(dBB,NOR,1);
for e = 1:numel(along)
    plot3(tC(:,1),tC(:,2),tC(:,3),'k.','MarkerSize',1)
    hold on
    fidx = find((lav < along(e) + dTH ) & (lav > along(e) - dTH));
    plot3(tC(fidx,1),tC(fidx,2),tC(fidx,3),'r.','MarkerSize',1)    
    subSet = tC(fidx,:);
    alongSub = PCA_REPROJ(subSet,dBB(:,e),BB(e,:));
    subSet = subSet - (dBB(:,e)*alongSub')';
    [aS aC aU(e,:) aE(:,:,e) gL gERR gLAM] = PCA_FIT_FULL(subSet,2);
    if e ~= 1
        for col = 1:size(E,2)
            if sign(aE(:,col,e-1)'*aE(:,col,e)) < 0
                aE(:,col,e) = -aE(:,col,e);
            end
        end
        
    end
    quiver3(BB(e,1),BB(e,2),BB(e,3),aE(1,1,e),aE(2,1,e),aE(3,1,e),4,'Color','r');
    quiver3(BB(e,1),BB(e,2),BB(e,3),aE(1,2,e),aE(2,2,e),aE(3,2,e),4,'Color','r');
    %quiver3(BB(e,1),BB(e,2),BB(e,3),NOR(1,e),NOR(2,e),NOR(3,e),4,'Color','b');
    %quiver3(BB(e,1),BB(e,2),BB(e,3),C(1,e),C(2,e),C(3,e),4,'Color','g');
    %hold off
    drawnow
end
plot3(BB(:,1),BB(:,2),BB(:,3),'g','LineWidth',3);
% project into frame for each point
mC = [];
for e = 1:size(tC,1)
    mC(e,:) = [PCA_REPROJ(tC(e,:),aE(:,:,aidx(e)),aU(aidx(e),:)) lav(e)];
end
%% see how basis changes with init tip angle
close all
for e = 1:size(aE,3)
    alongU = PCA_BKPROJ(BB(e,:),tE,tU);
    basisU1 = PCA_BKPROJ(BB(e,:) + squeeze(aE(:,1,e))',tE,tU);
    basisU2 = PCA_BKPROJ(BB(e,:) + squeeze(aE(:,2,e))',tE,tU);
    deltaU1 = basisU1 - alongU;
    deltaU2 = basisU2 - alongU;
    
    hold on
    plot(alongU,'r');
    plot(deltaU1,'g');
    plot(deltaU2,'b');
    drawnow
    pause(.1)
end
%% regress andPCA init T
close all
for tm =1:size(tipData,2)
    plot(tipData(:,1),tipData(:,tm),'.')
    drawnow
    pause(.2)
end
%%
close all
for e = 1:size(tipData,2)
    p(e,:) = polyfit(tipData(:,1),tipData(:,e),1);
    initF(:,e) = polyval(p(e,:),tipData(:,1));
end
delta = tipData - initF;
plot(delta');
[dS dC dU dE dL dERR dLAM] = PCA_FIT_FULL(delta,2);
figure;
plot(dE)
SIM = initF + dS;
figure;
for e = 1:50%size(tipData,1)
    plot(tipData(e,:),'k')
    hold on
    plot(SIM(e,:),'r')
    hold off
    axis([0 61 -pi/4 pi/2])
    drawnow    
    pause(.1)
end
uC = mean(dC,1);

dim = 2;
K = linspace(min(dC(:,dim)),max(dC(:,dim)));
figure;

for e = 1:numel(K)
    for tm = 1:61
        baseline(tm) = polyval(p(tm,:),mean(tipData(:,1)));
    end
    uC(dim) = K(e);
    modelF = PCA_BKPROJ(uC,dE,dU);
    modelF = modelF + baseline;
    plot(modelF);
    hold on
end
figure;plot(dE,'r--')
%% double regression
close all
DA = tipData(:,end) - tipData(:,1);
IA = tipData(:,1);
for e = 1:size(tipData,2)
    [vw(1) vw(2)] = view;
    subData = [DA,IA,tipData(:,e)];
    %[iS iC iU iE iL iERR iLAM] = PCA_FIT_FULL(subData,2);
    plot3(DA,IA,tipData(:,e),'.');
    %[icasig, A, W] = fastica(iC');
    %W = PCA_BKPROJ(W,iE,[0 0 0]);
    view(vw);
    drawnow
    pause(2)
end
%%
for tr = 1:size(toFit,1)
    myR(tr,:) = linspace(toFit(tr,1),toFit(tr,end),61);
end
%%
for e = 1:size(toFit,2)
    p(e,:) = polyfit(toFit(:,1),toFit(:,e),1);
    initF(:,e) = polyval(p(e,:),toFit(:,1));
end
close all
DEL = toFit - myR;
%[dS2 dC2 dU2 dE2 dL2 dERR2 dLAM2] = PCA_FIT_FULL(DEL,3);
kidx = kmeans(max(DEL,[],2),3);
nFit = bsxfun(@minus,toFit,toFit(:,1));
%nFit = bsxfun(@times,toFit,toFit(:,end).^-1);
CL = {'r.' 'g.' 'b.' 'k.' 'c.'};

VEC1 = [[0 pi/2];[pi/2 0]];
dVEC1 = diff(VEC1,1,2);
dVEC1 = dVEC1/norm(dVEC1);
dVEC1 = dVEC1*pi/2;
dVEC1 = VEC1(:,1) + dVEC1;
NCP = dVEC1;
LIN2 = [LIN(:,1) NCP];


LIN = [[-.2761 0];[-.2761 pi/2]];
CP = diff(LIN,1,2)/2;
LEN = norm(CP)
BI = CP / LEN;
CP = CP + LIN(:,1);
BI = [BI(2);-BI(1)];
BI = [CP CP + 100*BI];
%NCP = [1.06;.508];

RAD = LIN(:,2) - NCP;
OT = [LIN(:,2) NCP];
RAD = norm(RAD);
CTT = [];
for l = 1:100
    for tm = 1:size(toFit,2)
        
        plot(toFit(:,1),toFit(:,tm),'.')
        [iS iC iU iE iL iERR iLAM] = PCA_FIT_FULL([toFit(:,1),toFit(:,tm)],2);
        UFR{tm} = iU;
        EFR{tm} = iE;%*diag(diag(iLAM).^.5);
        LAMFR{tm} =iLAM;
        iC = iC*diag(diag(iLAM).^-.5);
        CTT = cat(3,CTT,iC);
        
        [iS2 iC2 iU2 iE2 iL2 iERR2 iLAM2] = PCA_FIT_FULL([toFit(:,end),toFit(:,tm)],2);
        iLAM2(iLAM2<0) = 0;
        UFR2{tm} = iU2;
        EFR2{tm} = iE2;%*diag(diag(iLAM).^.5);
        if tm > 1
            if sign(EFR2{tm-1}(:,1)'*EFR2{tm}(:,1)) < 1
                EFR2{tm}(:,1) = -EFR2{tm}(:,1);
            end
            if sign(EFR2{tm-1}(:,2)'*EFR2{tm}(:,2)) < 1
                EFR2{tm}(:,2) = -EFR2{tm}(:,2);
            end
        end
        LAMFR2{tm} =iLAM2;
        
        
        
        
        CIRC = [cos(TH);sin(TH)];
        CIRC = diag(diag(iLAM).^.5)*CIRC;
        CIRC = PCA_BKPROJ(CIRC',iE,iU);
        hold on
        plot(CIRC(:,1),CIRC(:,2),'r');
        %{
        for u = 1:5
            fidx = find(kidx==u);
            plot(toFit(fidx,1),nFit(fidx,tm),CL{u})
            hold on
        end
        %}
        
        Xi = linspace(min(toFit(:,1)),max(toFit(:,1)));
        L = polyval(p(tm,:),Xi);
        plot(Xi,L,'r')
        MCP = polyval(p(tm,:),0);
        
        val = ((pi/2)^2 + (pi/2)^2).^.5 - pi/2;
        TH = linspace(-pi,pi,100);
        CIR = (val + pi/2)*[cos(TH);sin(TH)];
        CIR2 = norm(RAD)*[cos(TH);sin(TH)];
        CIR2 = bsxfun(@plus,CIR2,NCP);
        CIR1 = bsxfun(@plus,CIR,[1;1]);
        CIR = bsxfun(@plus,CIR,[pi/2;0]);
        CIRP = pi/2*[cos(TH);sin(TH)];  
        CIRP1 = bsxfun(@plus,CIRP,dVEC1);
        CIRP = bsxfun(@plus,CIRP,[0;pi/2]);
        plot(0,MCP(1),'c*')
        plot(LIN(1,:),LIN(2,:),'c')
        plot(LIN2(1,:),LIN2(2,:),'r','LineWidth',4)
        plot(OT(1,:),OT(2,:),'r','LineWidth',4)        
        plot(CP(1),CP(2),'k*');
        plot(BI(1,:),BI(2,:),'c')
        plot(NCP(1),NCP(2),'*c')
        
        %plot(myR(:,1),myR(:,tm),'r.')
        plot(Xi,Xi,'g')
        plot([0 pi/2],[pi/2 0],'g')
        
        plot([-val pi/2],[0 0],'g')
        plot(CIR(1,:),CIR(2,:),'k')
        plot(CIR1(1,:),CIR1(2,:),'k')
        plot(CIR2(1,:),CIR2(2,:),'m*')
        plot(iU(1),iU(2),'y*')
        plot(CIRP(1,:),CIRP(2,:),'g')
        plot(CIRP1(1,:),CIRP1(2,:),'b')
        quiver(iU(1),iU(2),iE(1,1),iE(2,1),'k')
        quiver(iU(1),iU(2),iE(1,2),iE(2,2),'k')
        hold off
        
        %{
        plot(toFit(:,1),DEL(:,tm),'.')        
        %}
        axis([-pi/2 pi/2 -pi/2 pi/2])
        title(tm);
        pause(.1);
        
        drawnow
    end
    waitforbuttonpress
end
%%
close all
[SIM1 SIM2 SIM3 SIM4 SIM5 SIM6] = getSIM(toFit(:,1),toFit(:,end),UFR,EFR,LAMFR);
[SIM1_2 SIM2_2 SIM3_2 SIM4_2 SIM5_2 SIM6_2] = getSIM2(toFit(:,1),toFit(:,end),UFR2,EFR2,LAMFR2);
LEG{1} = 'regress--initTip time vs tip @ init--conditions start [init,init]';
LEG{2} = 'regress--initTip time vs tip @ init--conditions final [init,final]';
LEG{3} = 'regress--initTip time vs init_tip @ init--conditions final [init,final]';
LEG{4} = 'regress--initTip rev_time vs init_tip @ init--conditions final [init,final]';
for tr = 1:10%size(SIM,1)
    plot(SIM1(tr,:),'r');
    hold on
    plot(SIM2(tr,:),'b');    
    %plot(SIM4(tr,:),'m');
    %plot(SIM5(tr,:),'c');
    %plot(SIM6(tr,:),'ko');
    
    
    plot(SIM1_2(tr,:),'b--');
    hold on
    plot(SIM2_2(tr,:),'r--');    
    %plot(SIM4_2(tr,:),'m--');
    %plot(SIM5_2(tr,:),'c--');
    %plot(SIM6_2(tr,:),'ko');
    
    
    plot(toFit(tr,:),'k')
    axis([0 61 -pi/4 pi/2]);
    hold off
    drawnow
    %pause(1);
    legend(LEG)
    waitforbuttonpress
end
DEL1 = toFit - SIM1;
DEL2 = toFit - SIM2;
DEL3 = toFit - SIM1_2;
DEL4 = toFit - SIM1_2;
DEL1 = sum(DEL1.^2,2);
DEL2 = sum(DEL2.^2,2);
DEL3 = sum(DEL3.^2,2);
DEL4 = sum(DEL4.^2,2);
%%
%toFit = fMA;
close all
S = getPCAdy(toFit,[1 size(toFit,2)]);
NP = 20;
[n1 n2] = ndgrid(linspace(min(toFit(:,1)),max(toFit(:,1)),NP),linspace(min(toFit(:,end)),max(toFit(:,end)),NP));
for tm = 1:size(n1,2)
    SIM = getPCAdyf(S,[n1(:,tm) n2(:,tm)]);
    %dCg = polyval(p,n1(:,tm));
    %M = PCA_BKPROJ(dCg,dE,dU);
    %plot(SIM'+M','r');
    plot(SIM','b');
    %axis([0 61 -3*pi/4 3*pi/4]);
    axis([0 size(toFit,2) -40 120]);
    drawnow
    pause(1)
end
%%
R =readtext('/home/nate/Downloads/20141218_Spalding_all_individulas_raw.txt','\t');
H = R(1,:);
BLOCK = R(2:end,1:5);
%R{1,end+1} = 'lat1';
%R{1,end+1} = 'lat2';
%%
close all
toFit = cell2mat(R(2:end,5:113));
rmidx = any(abs(diff(toFit,1,2) > 5) ,2) | abs(toFit(:,1)) > 50;
%rmidx = any(abs(diff(toFit,1,2) > 5) ,2);
toFit(rmidx,:) = [];
nBLOCK = BLOCK;
nBLOCK(rmidx,:) = [];
plot(toFit');

%%
close all
SIM = getPCAdyf(S,toFit(:,[1 size(toFit,2)]));
DELTA = toFit - SIM;
[dS dC dU dE dL dERR dLAM] = PCA_FIT_FULL(DELTA,2);
p = polyfit(toFit(:,1),dC(:,1),1);
dCg = polyval(p,toFit(:,1));
%M = PCA_BKPROJ(dCg,dE,dU);
SIM2 = SIM + dS;
%SIM3 = SIM + M;
[sdC,pidx] = sort(dC(:,1));

for tr = 1:1:size(SIM,1)
    fidx = pidx(tr);
    plot(toFit(fidx,:),'k');
    hold on
    plot(SIM(fidx,:),'g');
    plot(SIM2(fidx,:),'b--');
    %plot(SIM3(tr,:),'g--');
    %axis([0 61 -pi/4 3*pi/4]);
    axis([0 size(toFit,2) -40 120]);
    hold off
    drawnow
    pause(.1);
end
%% glue for wolfgang

for e = 1:size(dC)
    nBLOCK{e,6} = dC(e,1);
    nBLOCK{e,7} = dC(e,2);
end
snBLOCK = nBLOCK;
snBLOCK(2:end+1,:) = nBLOCK(1:end,:);
H{6} = 'lat1';
H{7} = 'lat2';
for h = 1:7
    snBLOCK{1,h} = H{h};
end
snBLOCK(:,4:5) = [];
cell2csv('/home/nate/Downloads/20141218_Spalding_all_individulas_LAT.csv', snBLOCK);

%%
close all
X1 = linspace(min(toFit(:,1)),max(toFit(:,1)));
X2 = linspace(min(delta1(:,end)),max(delta1(:,end)));
[n1 n2] = ndgrid(X1,X2);
for tm = 1:size(toFit,2)
    F1 = polyval(p(tm,:),n1);
    F2 = polyval(q(tm,:),n2);
    FT = F1 + F2;
    
    [vw(1) vw(2)] = view;
    plot3(toFit(:,1),delta1(:,end),toFit(:,tm),'.')
    hold on
    mesh(n1,n2,FT);
    view(vw);
    hold off
    drawnow
    pause(.5);
end
%% regress and PCA init T double regression
close all
initF = [];
initG = [];
toFit = fMA;
%toFit = bsxfun(@minus,toFit,toFit(:,1));
for e = 1:size(toFit,2)
    p(e,:) = polyfit(toFit(:,1),toFit(:,e),1);
    initF(:,e) = polyval(p(e,:),toFit(:,1));
end
delta1 = toFit - initF;
for e = 1:size(toFit,2)
    xi = toFit(:,end)-toFit(:,1);
    xi = delta1(:,end);
    q(e,:) = polyfit(xi,delta1(:,e),1);
    initG(:,e) = polyval(q(e,:),xi);
end
plot(initG')
delta2 = delta1 - initG;
plot(delta2');
%delta2 = bsxfun(@minus,delta2,mean(delta2,2));
[dS2 dC2 dU2 dE2 dL2 dERR2 dLAM2] = PCA_FIT_FULL(delta2,3);
%[icasig, A, W] = fastica(dC2');
%K = bsxfun(@plus,A*icasig,mean(dC2,1)');
%modelF = PCA_BKPROJ(A,dE2,dU2);
plot(dE2)

dim = 1;
K = linspace(min(dC2(:,dim)),max(dC2(:,dim)));
figure;
uC2 = mean(dC2,1);
for e = 1:numel(K)
    for tm = 1:61
        baseline(tm) = polyval(p(tm,:),mean(toFit(:,1)));
        baseline2(tm) = polyval(q(tm,:),mean(xi));
    end
    uC2(dim) = K(e);
    modelF = PCA_BKPROJ(uC2,dE2,dU2);
    modelF = modelF + baseline + baseline2;
    plot(modelF);
    hold on
end
figure;
SIM = dS2 + initG + initF;
for e = 1:10%size(SIM,1)
    plot(SIM(e,:));
    hold on;
    plot(toFit(e,:),'k')
    
    hold off
    
    axis([0 61 -pi/4 pi/2])
    drawnow
    pause(.2);
end
%% plot populations for double regression
%kidx = kmeans(grData,5);
close all
UQ = unique(POP);
%UQ = unique(kidx);
dC2 = icasig';
for u = 1:numel(UQ)
    fidx = find(strcmp(POP,UQ{u}));
    %fidx = find(kidx == UQ(u));
    plot3(dC2(fidx,1),dC2(fidx,2),dC2(fidx,3),'.');
    hold all
    
end


%%
close all
X = [toFit(:,1),toFit(:,end) - toFit(:,1)];

UQ = unique(plateName);
TH = linspace(-pi,pi,100);

hold on
rmidx = zeros(size(toFit,1),1);
subT = [];
ASS = [];
AS = [];
for u = 1:numel(UQ)
    %plot(X(:,1),X(:,2),'.')
    %hold on
    
    fidx = strcmp(plateName,UQ{u});
    if sum(fidx) < 10
        rm(fidx) = 1;
    else
        %plot(X(fidx,1),X(fidx,2),'r.')
        [S C U E L ERR LAM] = PCA_FIT_FULL(X(fidx,:),2);
        if E(1,2) < 0
            E(:,2) = -E(:,2);
        end
        AS = [AS;atan(E(2,1)/E(1,1))];
        nC = (LAM.^.5)*CIR;    
        nC = PCA_BKPROJ(nC',E,U);
        subT = [subT;toFit(fidx,:)];
        ASS = [ASS;AS(end)*ones(sum(fidx),1)];
        %{
        plot(nC(:,1),nC(:,2),'r')
        plot(U(1),U(2),'k*')
        title(180/pi*AS(end));
        hold off
        drawnow
        %}
        
    end
    %waitforbuttonpress
    u
end
%%
UQ(find(rmidx)) = [];
pl = plateName;
pl(find(rmidx)) = [];

%%
close all
[mv,gidx] = min(abs(AS));
gidx = find(abs(ASS) == mv);
plot(subT(gidx,:)')

figure;
[mv,gidx] = max(AS);
gidx = find(ASS == mv);
plot(subT(gidx,:)')

figure;
[mv,gidx] = min(AS);
gidx = find(ASS == mv);
plot(subT(gidx,:)')





%% regress with init Tip angle
close all
%nTIP = bsxfun(@minus,tipData,mean(tipData,1));
nTIP = otipData + pi/2;
nTIP = cat(3,nTIP,repmat(nTIP(:,1),[1 size(nTIP,2)]));

%DL = dot(nTIP,nTIP,3).^-.5;
%nTIP = bsxfun(@times,nTIP,DL);
TH = atan(nTIP(:,:,1).*nTIP(:,:,2).^-1);
%{
ROT = [[cos(pi/4) sin(pi/4)];[-sin(pi/4) cos(pi/4)]];
for i = 1:size(nTIP,1)
    for j = 1:size(nTIP,2)
        nTIP(i,j,:) = ROT*squeeze(nTIP(i,j,:));
    end
    i
end
bsxfun(@times,nTIP,ROT)
nTIP = nTIP(:,:,1).*nTIP(:,:,2).^-2;
%nTIP = tipData;%-2*min(tipData(:));
%}
kidx = kmeans(TH,7);
CL = {'r.' 'g.' 'k.' 'b.' 'y.'};

for e = 2:size(tipData,2)
    
    [beta{e},PSI,stats,B{e}] = nlmefit(tipData(:,1),tipData(:,e),kidx,[],@(PHI,X)myRE(PHI,X),[1 0]);    
    for g = 1:5
        fidx = kidx == g;
        plot3(tipData(fidx,1),tipData(fidx,e),g*ones(sum(fidx),1),CL{g});
        
        axis([-.8 2 -2 2]);
        hold on
        
        
        hold on
        ti = linspace(min(tipData(:,1)),max(tipData(:,1)),100);
        tB = beta{e} + B{e}(:,g);
        yi = myRE(tB,ti);
        plot(ti,yi,CL{g}(1));
        
        
    end    
    hold off
    beta 
    drawnow
    %pause()
end
%%






%%
osubGravi = bsxfun(@minus,subGravi,subGravi(:,1));
%% backup data
tipDataBK = tipData;
specDataBK = specData;
predictionBK = prediction;
genoTypeBK = genoType;
plBK = pl;
%%
fidx = find(~strcmp(genoType,'null'));
fidx = [];
for e = 1:numel(genoType) 
    if ~strcmp(genoType{e}(1),'Z')
        if ~ strcmp(genoType{e},'null')
            fidx = [fidx e];
        end
        
    end
end

%% build and add starch
sd = readtext('/home/nate/Downloads/Starch-NIRdata.csv');
%% add starch
sY = cell2mat(sd(2:361,3:4));
sS = cell2mat(sd(2:361,6:end));
sS(:,end) = [];
[XL,YL,XS,YS,beta,PCTVAR] = plsregress(sS,sY,15);
sP = [ones(size(sS,1),1) sS]*beta;
prediction = [prediction [ones(size(specData,1),1) specData]*beta];
preTest = [preTest [ones(size(specData,1),1) specData]*beta];
%% add air space
prediction = [prediction prediction(:,end)-prediction(:,end-1)];
%% remove air space
prediction(:,end) = [];
%% remove starch
prediction(:,end-1:end) = [];
%%
plateNameBK = plateName;
%% correct WI08
for e = 1:numel(plateName)
    fidx = strfind(plateName{e},'WI08');
    if ~isempty(fidx)
        plateName{e} = [plateName{e}(1:4) '-' plateName{e}(5:end)];
        %[plateName{e}(1:4) '-' plateName{e}(5:end)]
    end
end
%% correct WI
for e = 1:numel(plateName)
    fidx = strfind(plateName{e},'WISN10');    
    if ~isempty(fidx)
        %plateName{e} = [plateName{e}(1:4) '-' plateName{e}(5:end)];
        plateName{e} = strrep(plateName{e},'-','_');
    end
end
%% correct WI11
for e = 1:numel(plateName)
    fidx = strfind(plateName{e},'WISN11');    
    if ~isempty(fidx)
        %plateName{e} = [plateName{e}(1:4) '-' plateName{e}(5:end)];
        plateName{e} = strrep(plateName{e},'-','_');
    end
end
%% correct WI09
for e = 1:numel(plateName)
    fidx = strfind(plateName{e},'WI09');
    if ~isempty(fidx)
        plateName{e} = [plateName{e}(1:4) '-' plateName{e}(5:end)];
        %[plateName{e}(1:4) '-' plateName{e}(5:end)]
    end
end
%% correct WI09
for e = 1:numel(plateName)
    plateName{e} = strrep(plateName{e},'-Nor','_Nor');
    plateName{e} = strrep(plateName{e},'-Mut','_Mut');
end
%% loook at pairs
for tr = 1:size(oldTip,1)
    plot(oldTip(tr,:),'b')
    hold on;plot(tipData(tr,:),'r')
    hold off 
    waitforbuttonpress
end
%% get lat para for NAM gravitropism
D = bsxfun(@minus,tipData,tipData(:,1));
%D = f.tipAngle;
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL(D,3);
[S1 C1 U1 E1 L1 ERR1 LAM1] = PCA_FIT_FULL(gradient(D),1);
MODEL = C1\wC;
SIM = C1*MODEL;
CS = linspace(min(C1),max(C1),10);
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM,wE,wU);
figure;
subplot(1,2,1);
plot(full_SIM');
DELTA = wC - SIM;
%DELTA = PCA_BKPROJ(DELTA,wE,wU);

[S2 C2 U2 E2 L2 ERR2 LAM2] = PCA_FIT_FULL(DELTA,1);
%check lat para
MODEL2 = C2\wC;
CS = linspace(min(C2),max(C2),10);
SIM2 = CS'*MODEL2;
SIM2 = PCA_BKPROJ(SIM2,wE,wU);
subplot(1,2,2);
plot(SIM2');
LAT = [C1 C2];
%% look at different populations
close all
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
UQ = unique(POP);
for u = 1:numel(UQ)
    fidx = strcmp(POP,UQ{u});
    otipData = bsxfun(@minus,tipData(fidx,:),tipData(fidx,1));    
    sU = mean(specData(fidx,:),1);
    sS = std(specData(fidx,:),1,1);%*sum(fidx).^-.5;
    
    
    pS = std(prediction(fidx,:));
    pU = mean(prediction(fidx,:));
    pS = pS.*pU.^-1
    
    
    tU = mean(otipData,1);
    tS = std(otipData,1,1);%*sum(fidx).^-.5;
    
    figure(h1)
    errorbar(tU*180/pi,tS*180/pi);
    hold all    
    
    figure(h2)
    errorbar(sU,sS);
    hold all    
    
    figure(h3);
    plot(tS*180/pi);
    hold all
    
    figure(h4);
    plot(pS);
    hold all
    NMS = {'weight','protein','oil','total_density','bone_density','bone_volume','total_volume','Starch_per','Starch_mg'};
    set(gca,'XTickLabel',NMS,'XTick',1:numel(NMS));
    
    hold all
    
end
figure(h1)
legend(UQ)
figure(h2)
legend(UQ)
figure(h3)
legend(UQ)
figure(h4)
legend(UQ)
%% cca
close all
fidx = 1:size(tipData,1);
otipData = bsxfun(@minus,tipData(fidx,:),tipData(fidx,1));
[utip UQ] =  generateUQmeans(otipData,pl);
[uspecData UQ] = generateUQmeans(specData,pl);
[uP UQ] =  generateUQmeans(prediction,pl);
RMSEP = [];
RMSEC = [];
RMSEC_s = [];
RMSEP_s = [];
for ns = 5:100
    [CORR masterX masterY mU mV RMSEPi RMSECi] = myHOM(uspecData,utip,pl,ns,5);
    RMSEP = [RMSEP mean(RMSEPi)];
    RMSEP_s = [RMSEP_s std(RMSEPi)];
    RMSEC = [RMSEC mean(RMSECi)];
    RMSEC_s = [RMSEC_s std(RMSECi)];
    errorbar(RMSEP,RMSEP_s,'b');
    plot(RMSEP);
    hold on
    errorbar(RMSEC,RMSEC_s,'r');
    hold off
    drawnow
    
end
%% cca
close all
fidx = 1:size(tipData,1);
otipData = bsxfun(@minus,tipData(fidx,:),tipData(fidx,1));
[utip UQ] =  generateUQmeans(otipData,pl);
[uspecData UQ] = generateUQmeans(specData,pl);
[uP UQ] =  generateUQmeans(prediction,pl);
PRESS = [];
for e = 1:size(specData,1)
    UQ{e} = num2str(e);
end
for ns = 1:130
    [PRESSi] = myHOMpls(specData,otipData,UQ,ns);    
    options = statset('UseParallel','always');
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(specData,otipData,ns,'Options',options,'cv',10);
    PRESS =  [PRESS mean(MSE,2)];
    plot(PRESS');
    drawnow
end
%% %% hold out via matlab commands tree
MRSE_train= [];
MRSE_test = [];
MRSE_test_itr = [];
MRSE_train_itr = [];
PVAL = [];
close all
q=1
X = specData;
Y = otipData;
%X = bsxfun(@minus,X,mean(X,1));
%Y = bsxfun(@minus,Y,mean(Y,1));
%[tS X tU tE tL tERR tLAM] = PCA_FIT_FULL(X,60);
%[sS Y sU sE sL sERR sLAM] = PCA_FIT_FULL(Y,60);

h1 = figure;
h2 = figure;
options = statset('UseParallel','always');
  

for factors = 1:100
    X = specData;
    Y = otipData;
    [tS Y tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);
    [sS X sU sE sL sERR sLAM] = PCA_FIT_FULL(X,30);
    
    func = @(xtrain,ytrain,xtest,ytest)myHoldOutTree(factors,xtrain,ytrain,xtest,ytest);
    tmpTrain = {};
    tmpTest = {};
    parfor itr = 1:10
        tic
        tmp = crossval(func,X,Y,'holdout',.50);
        tmp = reshape(tmp,[size(X,1) 2]);
        fidx0 = find(tmp(:,1)==0);
        fidx1 = find(tmp(:,1)==1);
        MRSE_train_itr(itr) = mean(tmp(fidx0,2));
        MRSE_test_itr(itr) = mean(tmp(fidx1,2));
        tmpTrain{itr} = tmp(fidx0,2)';
        tmpTest{itr} = tmp(fidx1,2)';
        [h(itr) p(itr)] = ttest2(tmp(fidx0,2),tmp(fidx1,2));
        toc
    end
    MRSE_train(:,factors) = [mean(MRSE_train_itr(:));std(MRSE_train_itr(:))];%.*numel(MRSE_train_itr).^-.5];
    MRSE_test(:,factors) = [mean(MRSE_test_itr(:));std(MRSE_test_itr(:))];%.*numel(MRSE_test_itr).^-.5];
    for itr = 1:size(tmpTrain)
        tmpTrain{itr} = mean(tmpTrain{itr});
        tmpTest{itr} = mean(tmpTest{itr});
    end
    [h PVAL(factors)] = ttest2(cell2mat(tmpTrain),cell2mat(tmpTest),[],[],'unequal');
    %PVAL(:,factors) = [mean(p);std(p)];
    figure(h1);
    errorbar(MRSE_train(1,:),MRSE_train(2,:),'r')
    hold on
    errorbar(MRSE_test(1,:),MRSE_test(2,:),'b')
    %errorbar(fG(1,:),fG(2,:),'g')
    
    
    drawnow
    figure(h2)
    %errorbar(PVAL(1,:),PVAL(2,:));
    plot(-log10(PVAL));
    hold on
    plot(-log10(.05*ones(size(PVAL))),'r');
    plot(-log10(.01*ones(size(PVAL))),'r');
    drawnow
end
%%
Ypre = [];
for e = 1:size(Y,2)
    ens = fitensemble(X,Y(:,e),'Bag',5,'Tree','type','regression');
    Ypre = [Ypre ens.predict(X)];
end
M = PCA_BKPROJ(Ypre,tE,tU);  
close all
plot(diag(corr(M,otipData)))
%%
ens = fitensemble(X,Y(:,1),'Bag',1,'Tree','type','regression','holdout',0.2);
plot(kfoldLoss(ens,'mode','cumulative'));
%% hold out via matlab commands GOOD HOLD OUT
MRSE_train= [];
MRSE_test = [];
MRSE_test_itr = [];
MRSE_train_itr = [];
PVAL = [];
close all
q=1
X = specData;
Y = otipData;
%X = bsxfun(@minus,X,mean(X,1));
%Y = bsxfun(@minus,Y,mean(Y,1));
%[tS X tU tE tL tERR tLAM] = PCA_FIT_FULL(X,60);
%[sS Y sU sE sL sERR sLAM] = PCA_FIT_FULL(Y,60);
X = uspecData;
Y = utip;
X = specData;
Y = otipData;
X = mU;
%[tS Y tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);
%[sS X sU sE sL sERR sLAM] = PCA_FIT_FULL(X,15);
h1 = figure;
h2 = figure;
options = statset('UseParallel','always');
  

for factors = 1:3
    func = @(xtrain,ytrain,xtest,ytest)myFUN(factors,xtrain,ytrain,xtest,ytest);
    tmpTrain = {};
    tmpTest = {};
    parfor itr = 1:100
        tic
        tmp = crossval(func,X,Y,'holdout',.20);
        tmp = reshape(tmp,[size(X,1) 2]);
        fidx0 = find(tmp(:,1)==0);
        fidx1 = find(tmp(:,1)==1);
        MRSE_train_itr(itr) = mean(tmp(fidx0,2));
        MRSE_test_itr(itr) = mean(tmp(fidx1,2));
        tmpTrain{itr} = tmp(fidx0,2)';
        tmpTest{itr} = tmp(fidx1,2)';
        [h(itr) p(itr)] = ttest2(tmp(fidx0,2),tmp(fidx1,2));
        toc
    end
    MRSE_train(:,factors) = [mean(MRSE_train_itr(:));std(MRSE_train_itr(:))];%.*numel(MRSE_train_itr).^-.5];
    MRSE_test(:,factors) = [mean(MRSE_test_itr(:));std(MRSE_test_itr(:))];%.*numel(MRSE_test_itr).^-.5];
    for itr = 1:size(tmpTrain)
        tmpTrain{itr} = mean(tmpTrain{itr});
        tmpTest{itr} = mean(tmpTest{itr});
    end
    [h PVAL(factors)] = ttest2(cell2mat(tmpTrain),cell2mat(tmpTest),[],[],'unequal');
    %PVAL(:,factors) = [mean(p);std(p)];
    figure(h1);
    errorbar(MRSE_train(1,:),MRSE_train(2,:),'r')
    hold on
    errorbar(MRSE_test(1,:),MRSE_test(2,:),'b')
    %errorbar(fG(1,:),fG(2,:),'g')
    
    
    drawnow
    figure(h2)
    %errorbar(PVAL(1,:),PVAL(2,:));
    plot(-log10(PVAL));
    hold on
    plot(-log10(.05*ones(size(PVAL))),'r');
    plot(-log10(.01*ones(size(PVAL))),'r');
    drawnow
end





%% genotype hold out X number of groups YES THIS ONE
close all

X = specData;
%X = mU;
%Y = tipData;
Y = otipData;
%Y = grData(:,end);


[S Y U Ex L ERR LAM] = PCA_FIT_FULL(Y,3);


%Y = bsxfun(@minus,Y,mean(Y,1));
%X = bsxfun(@minus,X,mean(X,1));
%{
for e = 1:size(Y,1)
    Y(e,:) = Y(e,:) * norm(Y(e,:))^-1;
end
for e = 1:size(X,1)
    X(e,:) = X(e,:) * norm(X(e,:))^-1;
end
%}

fG = [];
fG2 = [];
fG3 = [];

MASFG = [];
MASFG2 = [];
MAD = [];
UMT = [];
UMB = [];
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
h5 = figure;
for pull = 1:50
    fG = [];
    fG2 = [];
    fG3 = [];
    R = randi(numel(UQ),1,1);
    sidx0 = [];
    sidx1 = [];
    for e = 1:numel(R)
        %sidx0 = [sidx0 find(~strcmp(pl,UQ{R(e)}))];
        sidx1 = [sidx1 find(strcmp(pl,UQ{R(e)}))];
    end
    sidx0 = setdiff(1:numel(pl),sidx1);
    trainX = X(sidx0,:);
    trainY = Y(sidx0,:);
    testX = X(sidx1,:);        
    %testX = mean(testX,1);
    testY = Y(sidx1,:);
    %testY = mean(testY,1);

    uS = mean(trainX);
    [A,B,r,trainX,trainY,stats] = canoncorr(trainX,trainY);
    testX = bsxfun(@minus,testX,uS);
    testX = testX*A;
    
    options = statset('UseParallel','always');
    for factors = 3
        GENO_TEST = [];
        tic        
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(trainX,trainY,factors,'Options',options);
        gTest = [ones(size(testX,1),1),testX]*BETA;
        
        gTest2 = diag(corr(PCA_BKPROJ(gTest,Ex,U),PCA_BKPROJ(testY,Ex,U)));
        figure(h5);
        plot(gTest2);
        
        
%        [J sortIDX] = sort(gTest(:,50));
        [J sortIDX] = sort(gTest(:,end));
        
        PERCENT = .05;
        N = round(size(gTest,1)*PERCENT);
        UT = mean(testY(sortIDX(1:N),:),1);
        ST = std(testY(sortIDX(1:N),:),1,1);
        
        UB = mean(testY(sortIDX(end-N:end),:),1);
        SB = std(testY(sortIDX(end-N:end),:),1,1);
        
        UMT(factors,:,pull) = UT;
        UMB(factors,:,pull) = UB;
        
        
        
        
        gTest = gTest - testY;
        %gTest = sum(gTest.*gTest,2).^.5;        
        gTest = sum(gTest.*gTest.^2,2).^.5;
        %GENO_TEST= [GENO_TEST mean(gTest)];
        GENO_TEST = gTest;
        GENO_TEST2 = gTest2;
        toc

        fG = [fG [mean(GENO_TEST);std(GENO_TEST)*numel(GENO_TEST).^-.5]];
        fG2 = [fG2 [max(GENO_TEST2);std(GENO_TEST2)*numel(GENO_TEST2).^-.5]];
        %fG2 = [fG2 [mean(GENO_TEST2);std(GENO_TEST2)*numel(GENO_TEST2).^-.5]];
        fG3 = [fG3 GENO_TEST2];
        
        
        figure(h1);
        errorbar(fG(1,:),fG(2,:),'g')
        hold on
        errorbar(mean(MASFG,1),std(MASFG,1,1)*size(MASFG,1).^-.5,'r')        
        %errorbar(fG2(1,:),fG2(2,:),'r')
        hold off
        drawnow
        
        
        figure(h2);
        errorbar(fG2(1,:),fG2(2,:),'r')
        hold on;        
        %errorbar(mean(MASFG2,1),std(MASFG2,1,1)*size(MASFG2,1).^-.5,'g')
        errorbar(mean(MASFG2,1),std(MASFG2,1,1)*size(MASFG2,1).^-.5,'g')
        %errorbar(fG2(1,:),fG2(2,:),'r')
        hold off
        drawnow
        %{
        if size(MAD,3) > 1
            figure(h3);
            mesh(mean(MAD,3));
        end
        %}
        
        
    end
    
    
    figure(h4);
    N = round(size(otipData,1)*PERCENT);
    [J sIDX] = sort(otipData(:,50));
    UTOP = otipData(sIDX(1:N),:);
    UBOT = otipData(sIDX(end-N:end),:);
    plot(180/pi*mean(UTOP,1),'r')
    hold on
    plot(180/pi*mean(UBOT,1),'b')
    plot(180/pi*mean(UMT(end,:,:),3),'m');
    plot(180/pi*mean(UMB(end,:,:),3),'c');
    hold off
    
    MASFG = [MASFG;fG(1,:)];
    MAD = cat(3,MAD,fG3);
    MASFG2 = [MASFG2;fG2(1,:)];
end
%%
close all
figure;
errorbar(mean(MASFG,1),std(MASFG,1,1)*size(MASFG,1).^-.5,'r')
figure;
errorbar(mean(MASFG2,1),std(MASFG2,1,1)*size(MASFG2,1).^-.5,'g')
%% groups - GROUPS
%{
X = otipData;
NG = 3;
kidx = kmeans(X,NG);
kUQ = unique(kidx);
for u = 1:numel(kUQ)
    errorbar(mean(X(kidx==u,:),1),std(X(kidx==u,:)))
    hold all
end
waitforbuttonpress
close all
%}
close all

otipData = bsxfun(@minus,tipData,tipData(:,1));
X = specData;
X = mU;
%X = prediction;
Y = otipData;
%Y = grData;
%X = prediction;

%{
uX = [];
uY = [];
UQ = unique(pl);
for u = 1:numel(UQ)
    fidx = find(strcmp(pl,UQ{u}));
    uX(u,:) = mean(X(fidx,:));
    uY(u,:) = mean(Y(fidx,:));
end
X = uX;
Y = uY;
%}
NG = 3;
kidx = kmeans(Y(:,end),NG);
%kidx = kmeans(grData(:,end),NG);
KL = [];
for u = 1:NG
    KL(u) = mean(sum(Y(kidx==u,end),2));
end
[KL sidx] = sort(KL,'descend');
for u = 1:NG
    fidx = find(kidx==u);
    kkidx(fidx) = find(sidx == u);
end
kidx = kkidx;
CMtest = [];
CMtrain = [];

h1 = figure;
h2 = figure;
h3 = figure;
% display groups
CL = {'r' 'g' 'b' 'm' 'c'};
UU = [];
SU = [];
for u = 1:NG
    fidx = find(kidx==u);
    figure(h1);
    errorbar(180/pi*mean(otipData(fidx,:)),180/pi*std(otipData(fidx,:),1,1)/40,CL{u});
    UU = [UU;mean(grData(fidx,:))];
    SU = [SU;std(grData(fidx,:),1,1)*numel(fidx)^-.5];
    hold all
    figure(h2);
    errorbar(9.64*mean(grData(fidx,:)),9.64*std(grData(fidx,:),1,1)*size(fidx,2)^-.5,CL{u});
    hold all
end
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/UU_GR.csv',[UU;SU]);


%[tS X tU tE tL tERR tLAM] = PCA_FIT_FULL(X,50);
%X = [X bsxfun(@minus,prediction,mean(prediction,1))];

%X = X*diag(diag(tLAM).^-.5);
%X = [X Y(:,end)];
[sS Y sU sE sL sERR sLAM] = PCA_FIT_FULL(Y,4);


SELEC = [];
SELEC2 = [];
UQ = unique(pl);
for pull = 1:500
    pull
    %P = randperm(size(X,1));
    %X = X(P,:);
    
    %{
    NT = 1000;
    R = randperm(size(X,1));
    sidx1 = R(1:NT);
    sidx0 = R(NT+1:end);
    %}
    
    
    R = randi(numel(UQ),1,50);
    sidx0 = [];
    sidx1 = [];
    for e = 1:numel(R)        
        sidx1 = [sidx1 find(strcmp(pl,UQ{R(e)}))];
    end    
    sidx0 = setdiff(1:numel(pl),sidx1);
    
    
    
    trainX = X(sidx0,:);
    trainY = sS(sidx0,:);
    gtrainY = kidx(sidx0);
    testX = X(sidx1,:);        
    %testX = mean(testX,1);
    testY = sS(sidx1,:);
    gtestY = kidx(sidx1);
    testGR = grData(sidx1,:);
    
    
    %mdl = ClassificationKNN.fit(trainX,gtrainY,'NumNeighbors',15,'distance','mahalanobis');
    %[testLabel2,score] = predict(mdl,testX);
    
    %nbGau = NaiveBayes.fit(trainX, gtrainY);%;,'Distribution','kernel','KSWidth',.10);        
    %testLabel2 = nbGau.predict(testX);
    
    %tree = ClassificationTree.fit(trainX,gtrainY);
    %testLabel2 = predict(tree,testX);
    
    %cl = ClassificationDiscriminant.fit(trainX,gtrainY);
    %testLabel2 = predict(cl,testX);
    %{
    tmp = {};
    for u = 1:numel(gtrainY)
        tmp{u} = num2str(gtrainY);
    end
    ens = classregtree(trainX,tmp);
    testLabel = eval(ens,testX);
    for e = 1:numel(testLabel)
        testLabel{e} = str2num(testLabel{e});
    end
    testLabel2 = cell2mat(testLabel);
    %}
    %{
    net = patternnet(10);
    G = [];
    for u = 1:NG
        fidx = find(gtrainY==u);
        G(fidx,u) = 1;
    end
    gtrainY = G;
    net = train(net,trainX',gtrainY');
    testLabel2 = net(testX');
    testLabel2 = vec2ind(testLabel2)';    
    gtrainY = vec2ind(gtrainY')';
    %}
    
    %ens = fitensemble(trainX,gtrainY,'Bag',100,'Discriminant','type','classification');
    %testLabel2 = predict(ens,testX);
    
    %{
    kidx1 = find(testLabel2 == 1);
    SELEC = [SELEC;mean(testY(kidx1,:),1)];
    
    kidx1 = find(testLabel2 == 3);
    SELEC2 = [SELEC2;mean(testY(kidx1,:),1)];
    %}
    
    
    
    
    kidx1 = find(testLabel2 == 1);
    SELEC = [SELEC;mean(testGR(kidx1,:),1)];
    
    kidx1 = find(testLabel2 == 3);
    SELEC2 = [SELEC2;mean(testGR(kidx1,:),1)];
    
    
    
    
    
    
    %SELEC = [SELEC;mean(grData(kidx1,:),1)];
    
    %testLabel2 = nbGau.predict(trainX);
    %testLabel3 = classify(testX,trainX,gtrainY,'mahalanobis');
    
    %cl = ClassificationDiscriminant.fit(trainX,gtrainY);
    %testLabel6 = predict(cl,testX);
    %ens = fitensemble(trainX,gtrainY,'Bag',100,'Discriminant','type','classification');
    %testLabel7 = predict(ens,testX);
    
    %SVMStruct = svmtrain(trainX,gtrainY,'kernel_function','rbf','rbf_sigma',500);
    %testLabel4 = svmclassify(SVMStruct,testX);
    
    %testLabel = mode([testLabel1 testLabel2 testLabel3],2);
    %testLabel = round(testLabel);
    
    
    testLabel = testLabel2;
    
    
    
    %mdl = ClassificationKNN.fit(trainX,gtrainY);
    %[testLabel,score] = predict(mdl,testX);
    %[trainLabel,score] = predict(mdl,trainX);
    %SVMStruct = svmtrain(trainX,gtrainY,'kernel_function','rbf','rbf_sigma',500);
    %testLabel = svmclassify(SVMStruct,testX);
    %nbGau = NaiveBayes.fit(trainX, gtrainY);
    %trainLabel = nbGau.predict(trainX);
    %testLabel = nbGau.predict(testX);
    %testLabel = classify(testX,trainX,gtrainY);
    %ens = fitensemble(trainX,gtrainY,'Bag',500,'tree','type','classification');
    %testLabel = predict(ens,testX);
    
    %{
    net = patternnet(10);
    G = [];
    for u = 1:NG
        fidx = find(gtrainY==u);
        G(fidx,u) = 1;
    end
    gtrainY = G;
    net = train(net,trainX',gtrainY');
    testLabel = net(testX');
    testLabel = vec2ind(testLabel)';    
    gtrainY = vec2ind(gtrainY')';
    %}
    
    
    
    
    
    %{
    ens = classregtree(trainX,num2str(gtrainY));
    testLabel = eval(ens,testX);
    for e = 1:numel(testLabel)
        testLabel{e} = str2num(testLabel{e});
    end
    testLabel = cell2mat(testLabel);
    %}
    currCM = confusionmat(gtestY,testLabel);
    CMtest = cat(3,CMtest,currCM'); %/ numel(gtestY));
    %CMtrain = cat(3,CMtrain,confusionmat(trainLabel,gtrainY));% / numel(gtrainY));
    mean(CMtest,3);
    %mean(CMtrain,3)
   
    CL = {'m' 'c' 'r' 'b' 'g'};
    for u = 1:NG
        fidx = find(testLabel==u);
        errorbar(mean(testY(fidx,:)),std(testY(fidx,:)).*size(testY,1)^-.5,CL{u});
        hold on
        fidx = find(gtrainY==u);
        errorbar(mean(trainY(fidx,:),1),std(trainY(fidx,:)).*size(trainY,1)^-.5,'k');
        
    end
    %waitforbuttonpress
    drawnow
    pause(.1)
    close all
     
    drawnow
    %waitforbuttonpress
    %close all
    
    %JK = sum(CMtest,1);    
    %uJ = bsxfun(@times,CMtest,JK.^-1);
    %PERCM = mean(uJ,3)
    %PERCM_S = std(uJ,1,3);
    uJ = mean(CMtest,3);    
    PERCM = bsxfun(@times,uJ,sum(uJ,1).^-1)
    
    
    %figure(h3);hold off;errorbar(180/pi*nanmean(SELEC,1),180/pi*nanstd(SELEC,1,1)*size(SELEC,1)^-.5);hold all;errorbar(180/pi*nanmean(SELEC2,1),180/pi*nanstd(SELEC2,1,1)*size(SELEC2,1)^-.5);
    
    %figure(h3);hold off;errorbar(9.64*nanmean(SELEC,1),9.64*nanstd(SELEC,1,1)*size(SELEC,1)^-.5);hold all;errorbar(9.64*nanmean(SELEC2,1),9.64*nanstd(SELEC2,1,1)*size(SELEC2,1)^-.5);
    
end
%[p h] = ttest2(SELEC,SELEC2,[],[],[],1);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/prediction_growthRate.csv',[nanmean(SELEC,1);nanstd(SELEC,1,1);nanmean(SELEC2,1);nanstd(SELEC2,1,1)]);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/CM.csv',PERCM);
%% growth rate correlation to oil, starch etc
X = specData;
%X = prediction;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,grData(:,end),15);
SIM = [ones(size(X,1),1),X]*BETA;
corr(SIM,grData(:,end))
%%
nbGau = NaiveBayes.fit(X,kidx);%;,'Distribution','kernel','KSWidth',.10);    
testLabel2 = nbGau.predict(X);
h = confusionmat(kidx,testLabel2);
%% plot prediction data RAW
close all
%Y = gradient(otipData);
%mean(prediction(:,[2 4 5]),2).*mean(prediction(:,[1 6 7]),2).^-1
[raw_corr raw_pre] = corr([prediction prediction(:,1).*prediction(:,3) prediction(:,1).*prediction(:,2) prediction(:,1).^-1 -prediction(:,[2 4 5]) -prediction(:,[5]).*prediction(:,1).^-1],(otipData));figure;
figure;
plot(raw_corr'.*(raw_pre' < .01));
figure;
plot(raw_corr')
NMS = {'weight','protein','oil','total_density','bone_density','bone_volume','total_volume','Starch1','Starch2','Total_oil','Total_PRO','inverseW','POS1','POS2','POS3','RE'};
legend(NMS);
hold on;
CL = {'r' 'b' 'g' 'k' 'm' 'c' 'r' 'b' 'g' 'k' 'm'};
figure;
for e = 1:size(raw_corr,1)
    fidx = (raw_pre(e,:) < .01);
    R = regionprops(fidx,'PixelIdxList');
    for n = 1:numel(R)
        plot(R(n).PixelIdxList,raw_corr(e,R(n).PixelIdxList),CL{e});
        hold on;
    end
end
%% plot prediction data RAW 2
[raw_corr raw_pre] = corr(preTest,otipData);figure;
plot(raw_corr')


%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/raw_corr.csv',raw_corr);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/raw_pval.csv',raw_pre);
%%
close all
X = specData;
Y = otipData;
[tS X tU tE tL tERR tLAM] = PCA_FIT_FULL(X,5);
[sS Y sU sE sL sERR sLAM] = PCA_FIT_FULL(Y,3);
preY = [];
for e = 1:size(Y,2)
    R = classregtree(X,Y(:,e));
    preY = [preY eval(R,X)];
end
preY = PCA_BKPROJ(preY,sE,sU);
plot(diag(corr(preY,otipData)));
%% master plot
close all

plot(fG(1,:),'g');
hold on
plot(fG(1,:)+fG(2,:),'g--');
plot(fG(1,:)-fG(2,:),'g--');
plot(MRSE_train(1,:),'r');
plot(MRSE_train(1,:)+MRSE_train(2,:),'r--');
plot(MRSE_train(1,:)-MRSE_train(2,:),'r--');
plot(MRSE_test(1,:),'b');
plot(MRSE_test(1,:)+MRSE_test(2,:),'b--');
plot(MRSE_test(1,:)-MRSE_test(2,:),'b--');
csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/genoTypeHoldout.csv',fG);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/et_internalHoldout.csv',MRSE_train);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/et_externalHoldout.csv',MRSE_test);
%% genotype hold out single hold out
close all

X = specData;
Y = tipData;
fG = [];
fG2 = [];


for factors = 1:50
    R = randi(numel(UQ),1,50);
    GENO_TEST = [];
    GENO_TEST2 = [];
    parfor u = 1:numel(R)
        tic        
        sidx0 = find(~strcmp(pl,UQ{R(u)}));
        sidx1 = find(strcmp(pl,UQ{R(u)}));        
        trainX = X(sidx0,:);
        trainY = Y(sidx0,:);
        
        testX = X(sidx1,:);        
        %testX = mean(testX,1);
        testY = Y(sidx1,:);
        
        %testY = mean(testY,1);
        
        %{
        gTest = [];
        for e = 1:size(Y,2)
            ens = fitensemble(trainX,trainY(:,e),'Bag',factors,'Tree','type','regression');
            gTest = [gTest ens.predict(X)];
        end
        %}
        
    
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(trainX,trainY,factors);
        gTest = [ones(size(testX,1),1),testX]*BETA;
    
        gTest = gTest - testY;
        gTest2 = diag(corr(gTest,testY));
        gTest = sum(gTest.*gTest,2).^.5;
        GENO_TEST = [GENO_TEST mean(gTest)];
        GENO_TEST2 = [GENO_TEST2 mean(gTest2)];
        toc
    end
    fG = [fG [mean(GENO_TEST);std(GENO_TEST)*numel(GENO_TEST).^-.5]];
    fG2 = [fG2 [mean(GENO_TEST2);std(GENO_TEST2)*numel(GENO_TEST2).^-.5]];
    errorbar(fG(1,:),fG(2,:),'g')
    hold on
    errorbar(fG2(1,:),fG2(2,:),'r')
    hold off
    drawnow
end
%%
close all

X = specData;
Y = tipData;
fG = [];
fG2 = [];
R = randi(numel(UQ),1,50);
sidx0 = [];
sidx1 = [];
for e = 1:numel(R)
    sidx0 = [sidx0 find(~strcmp(pl,UQ{R(e)}))];
    sidx1 = [sidx1 find(strcmp(pl,UQ{R(e)}))];
    
end
trainX = X(sidx0,:);
trainY = Y(sidx0,:);
testX = X(sidx1,:);        
%testX = mean(testX,1);
testY = Y(sidx1,:);
%testY = mean(testY,1);

options = statset('UseParallel','always');
for factors = 1:50
    GENO_TEST = [];
    GENO_TEST2 = [];
    tic        
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(trainX,trainY,factors,'Options',options);
    gTest = [ones(size(testX,1),1),testX]*BETA;
    gTest = gTest - testY;
    gTest = sum(gTest.*gTest,2).^.5;
    %GENO_TEST= [GENO_TEST mean(gTest)];
    GENO_TEST = gTest;
    
    
    
    
    gTest2 = diag(corr(gTest,testY));
    GENO_TEST2 = [GENO_TEST2 mean(gTest2)];
    toc
    
    fG = [fG [mean(GENO_TEST);std(GENO_TEST)*numel(GENO_TEST).^-.5]];
    fG2 = [fG2 [mean(GENO_TEST2);std(GENO_TEST2)*numel(GENO_TEST2).^-.5]];
    errorbar(fG(1,:),fG(2,:),'g')
    hold on
    errorbar(fG2(1,:),fG2(2,:),'r')
    hold on
    %errorbar(fG2(1,:),fG2(2,:),'r')
    hold off
    drawnow
end

%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/MRSE_train.csv',MRSE_train);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/MRSE_test.csv',MRSE_test);

%%
X = specData;
Y = otipData;
[sS X sU sE sL sERR sLAM] = PCA_FIT_FULL(X,20);
M = crM(X,Y(:,end),{},1,2);
%% focused on corr hold out - THIS ONE TOO
close all
%[X,muX,sigmaX] = zscore(uspecData);
%[Y,muY,sigmaY] = zscore(utip);
X = uspecData;
Y = utip;
X = specData;
Y = otipData;
%X = prediction;

%[tS Y tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);
%[sS X sU sE sL sERR sLAM] = PCA_FIT_FULL(X,20);

U = [];
S = [];
cnt = 1;
figure;
for att = 50
   
    CT = [];
    PT = [];
    for e = 1%:50

        %R = randi(numel(UQ),1,621);
        R = randperm(numel(UQ));
        R = R(1:20);
        sidx0 = [];
        sidx1 = [];
        for i = 1:numel(R)
            sidx1 = [sidx1 find(strcmp(pl,UQ{R(i)}))];
        end
        %sidx0 = unique(sidx0);
        sidx0 = setdiff(1:numel(pl),sidx1);
        trainX = X(sidx0,:);
        trainY = Y(sidx0,:);
        testX = X(sidx1,:);
        testY = Y(sidx1,:);
        
        
        %{
        net = feedforwardnet([att 3 att]);
        net = train(net,trainX',trainY');
        Ypredict = net(testX')';
        %}
        
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(trainX,trainY,att);
        Ypredict = [ones(size(testX,1),1),testX]*BETA;
        
        
        
        %{
        Ypredict = [];
        for i = 1:size(trainY,2)
            ens = fitensemble(trainX,trainY(:,i),'Bag',att,'Tree','type','regression');
            %ens = RegressionTree.fit(trainX,trainY(:,i));
            %ens = crossval(ens);
            %ens = compact(ens);
            Ypredict = [Ypredict predict(ens,testX)];
            %Ypredict = [Ypredict ens.predict(testX)];
        end
        %}
        
        
        %Ypredict = PCA_BKPROJ(Ypredict,tE,tU);
        
        
                
        [C P] = corr(otipData(sidx1,:),Ypredict);
        C = diag(mean(abs((otipData(sidx1,:) - Ypredict)),1));
        
        CT(e,:) = diag(C)';
        PT(e,:) = diag(P)';
        
        for tr = 1:size(Ypredict,1)
            plot(180/pi*Ypredict(tr,:),'r')
            hold on
            plot(180/pi*(Ypredict(tr,:)+mean(CT,1)),'r--')
            plot(180/pi*(Ypredict(tr,:)-mean(CT,1)),'r--')
            plot(180/pi*otipData(sidx1(tr),:),'b')
            hold off
            axis([0 61 -10 90]);
            drawnow
            pause(.5)
        end
        
        
        
        e
    end
    U = [U;mean(CT,1)];
    S = [S;std(CT,1,1)*size(CT,1)^-.5];
    %errorbar(U(cnt,:),S(cnt,:));
    plot(U(cnt,:))
    cnt = cnt  +1;
    hold all
    drawnow
end
figure;
errorbar(U,S)
%hold on
%plot(mean(PT,1),'r')

%% make genotype error bar plot
figure;
for e = 1:size(U,1)
    errorbar(U(e,:),S(e,:))
    hold all
end
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/genoType_repeatHoldout_mean.csv',U);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/genoType_repeatHoldout_se.csv',S);
%%
figure;plot(Ypredict');
RMSE = (Ypredict-otipData);
    RMSE = sum(RMSE.*RMSE,2).^.5;
    title(num2str(mean(RMSE)));
    %{
    for e = 1:size(M,1)
        plot(otipData(e,:)*180/pi,'b');
        hold on
        plot(Ypredict(e,:)*180/pi,'r');
        plot(mean(otipData,1)*180/pi,'k');
        axis([0 61 -30 60]);
        drawnow
        hold off
        pause(.3)
    end
    %}
%%
X = uspecData;
Y = utip;
X = specData;
Y = otipData;
[tS Y tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);
[sS X sU sE sL sERR sLAM] = PCA_FIT_FULL(X,15);

%net = cascadeforwardnet([15 3 15]);
net = feedforwardnet([8 3 8]);
%net = feedforwardnet([15 8 15]);
%net = fitnet([3 8 3]);
net = train(net,X',Y');
preY = net(X');
M = PCA_BKPROJ(preY',tE,tU);
close all
plot(diag(corr(M,otipData)));
figure;plot(M');
RMSE = (M-otipData);
RMSE = sum(RMSE.*RMSE,2).^.5;
title(num2str(mean(RMSE)));
for e = 1:size(M,1)
        plot(otipData(e,:)*180/pi,'b');
        hold on
        plot(M(e,:)*180/pi,'r');
        plot(mean(otipData,1)*180/pi,'k');
        axis([0 61 -30 60]);
        drawnow
        hold off
        pause(.3)
    end
%%
%Ypredict = bsxfun(@times,Ypredict,sigmaY.^-1);
%Ypredict = bsxfun(@plus,Ypredict,muY);
mXCORR = corr(XS,uP);
mYCORR = corr(XS,utip);
[mPCORR pVAL] = corr(uP,utip);
for e = 1:size(XS,2)
    figure;
    bar(mXCORR(:,e));
    NMS = {'weight','protein','oil','total_density','bone_density','bone_volume','total_volume','Starch_per','Starch_mg'};
    set(gca,'XTickLabel',NMS,'XTick',1:numel(NMS));
end

figure;
for e = 1:size(XS,2)    
    plot(mYCORR(e,:));    
    hold all
    LEG{e} = num2str(e);
end
legend(LEG)


figure;
CL = {'r' 'b' 'g' 'k' 'c' 'm' 'y'  'r*' 'b*' 'g*' 'k*' 'c*' 'm*' 'y*'};
for e = 1:size(uP,2)    
    plot(mPCORR(e,:),CL{e});
    hold all
    LEG{e} = NMS{e};
end

for e = 1:size(uP,2)        
    plot(.05*(pVAL(e,:) < .01),CL{e});
    hold all
    LEG{e} = NMS{e};
end
legend(LEG)



%%
figure;
%plot(corr(.8*uP(:,2)-.4*uP(:,4),utip))
plot(corr(uP(:,2),utip))
%%
for e = 1:size(XS,2)
    figure;
    plot(XS(:,e),YS(:,e),'.');
    title(num2str(corr(XS(:,e),YS(:,e))));            
    hold on
    plot(linspace(-max(XS(:,e)),max(XS(:,e)),2),linspace(-max(XS(:,e)),max(XS(:,e)),2))
end
%%
close all
[CORR masterX masterY mU mV RMSEPi RMSECi] = myHOM(uspecData,utip,pl,5,5);
for e = 1:1%size(mU,2)
        figure;
        plot(mU(:,e),mV(:,e),'.');
        title(num2str(corr(mU(:,e),mV(:,e)))); 
        hold on
        plot(linspace(-4,4,2),linspace(-4,4,2))
    end
%%
close all
ccn = 2;
figure;
plot(corr(mV(:,ccn),utip)');
figure;
plot(corr(mU(:,ccn),utip)');
figure;
bar(corr(mU(:,ccn),uP));
NMS = {'weight','protein','oil','total_density','bone_density','bone_volume','total_volume','Starch_per','Starch_mg'};
set(gca,'XTickLabel',NMS,'XTick',1:numel(NMS));


%% analysis of spec to gravi via cca
DS = 1;

fidx = 1:size(tipData,1);
%fidx = strcmp(POP,'NAM_parents');
%fidx = strcmp(POP,'NAM_parents') | strcmp(POP,'NC-350_RILs');

%fidx = strcmp(POP,'NC-350_RILs');
%fidx = strcmp(POP,'seedling_phenotyping_widiv');
fidx = find(fidx);
%fidx = fidx(end-110:end);
%{
UQ = unique(pl(fidx));
%fidx=  [];

for u = 1:29
    f = find(strcmp(pl,UQ{u}));
    fidx = [fidx f];
end
%}
otipData = bsxfun(@minus,tipData(fidx,:),tipData(fidx,1));
%otipData = gradient(tipData);
%[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(otipData(:,1:DS:end),3);
%[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(specData(fidx,1:779),15);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(otipData,5);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(specData,15);


%{
gidx = strcmp(POP,'seedling_phenotyping_widiv');
gidx = strcmp(POP,'NC-350_RILs');
otipData = bsxfun(@minus,tipData(gidx,:),tipData(gidx,1));
tC = PCA_REPROJ(otipData(:,1:DS:end),tE,tU);
sC = PCA_REPROJ(specData(gidx,1:DS:end),sE,sU);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(otipData(:,1:DS:end),3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(specData(gidx,1:779),15);
mU = sC*mA;
mV = tC*mB;
%}
%[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(prediction,9);
close all
%X = [sC specData(:,780)-mean(specData(:,780))];
%X = [sC prediction(:,1)-mean(prediction(:,1))];
X = sC;
Y = tC;
K = sC\tC;

%X = X(randperm(size(X,1)),:);
%[X,muX,sigmaX] = zscore(X);

%[Y,muY,sigmaY] = zscore(Y);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
[mA,mB,mU,mV] = myCCA(X,Y,5);
for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    %set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Shape:' num2str(e)]);
end

[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mV,2)
    figure;
    bar(Ycore(e,:));
    %set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi));
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mU,2)
    figure;
    bar(Xcore(e,:));
    %set(gca,'XTickLabel',[LAB.spec],'XTick',1:numel([LAB.spec]))
    title(['Core Spec:' num2str(e)]);
end

for e = 1:size(mU,2)
    figure;
    plot(mU(:,e),mV(:,e),'.');
    [COR(e),R(e)] = corr(mU(:,e),mV(:,e));
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end


preY = (pinv(mB')*mU')';
%preY = sC*K;
for e = 1:size(Y,2)
    %[beta,SIGMA,RESID,COVB] = mvregress(X,Y(:,e));
    %beta = robustfit(X,Y(:,e));
    
    %[B,FitInfo] = lasso(X,Y(:,e),'Alpha',.9,'CV',10);
    
    %preY(:,e) = [ones(size(X,1),1) X]*beta;
    %preY(:,e) = X*B(:,FitInfo.IndexMinMSE);
end

%preY = bsxfun(@times,sigmaY.^-1,preY);
%p = polyfit(preY,Y,1);
%preY = polyval(p,preY);
corr(preY,Y);
M = PCA_BKPROJ(preY,tE,tU);

%{
delta = (otipData-M);
%delta = Y - preY;
[dS dC dU dE dL dERR dLAM] = PCA_FIT_FULL(delta,5);
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL(specData,10);
K = wC\dC;
[mA,mB,mr,mU,mV,mstats] = canoncorr(wC,dC);
preD = (inv(mB')*mU')';
preD = wC*K;
Md = PCA_BKPROJ(preD,dE,dU);
%G = preY+Md;
%G = PCA_BKPROJ(G,tE,tU);
G = M + Md;
%}
%delta = (otipData-(M+Md));
figure;
plot(corr(mV(:,1),otipData)');
figure;
plot(corr(mU(:,1),otipData)');
figure;
plot(corr(prediction(fidx,1),otipData)','r');
figure;
bar(corr(mV(:,1),prediction(fidx,1:end)));
NMS = {'weight','protein','oil','total_density','bone_density','bone_volume','total_volume','Starch_per','Starch_mg'};
set(gca,'XTickLabel',NMS,'XTick',1:numel(NMS));
figure;plot(diag(corr(M,otipData)));
%%
figure;
for e = 1:size(M,1)
    plot(180/pi*M(e,:),'r');
    hold on
    plot(180/pi*(G(e,:)),'c');
    hold on
    plot(180/pi*otipData(e,1:DS:end),'b');
    plot(180/pi*tS(e,:),'k');
    %plot(mean(tipData,1)*180/pi,'k')
    hold off
    axis([0 61 -20 90]);
    pause(.3)
end
%%
delta = otipData - M;
delta = sum(delta.*delta,2).^.5;
fidx = delta < 2;
figure;plot(mU(:,1),mV(:,1),'.')
hold on;plot(mU(fidx,:),mV(fidx,1),'r.');

figure;
fidx = find(fidx);
figure;plot(diag(corr(M(fidx,:),otipData(fidx,:))));
figure;
for i = 1:numel(fidx)
    e = fidx(i);
    plot(180/pi*M(e,:),'r');    
    hold on
    plot(180/pi*otipData(e,1:DS:end),'b');
    plot(180/pi*tS(e,:),'k');
    %plot(mean(tipData,1)*180/pi,'k')
    hold off
    axis([0 61 -20 90]);
    pause(.3)
end
%%
specDataBK = specData;
%% 
specData = specDataBK;
%%
close all
PRE = [];
%otipData = gradient(otipData);
specData = specData(fidx,:);
for e = 1:size(sC,1)
    %sub = specData(setdiff(1:size(sC,1),e),:);
    %subT = otipData(setdiff(1:size(sC,1),e),:);    
    delta = bsxfun(@minus,specData,specData(e,:));
    delta = sum(delta.*delta,2);
    %delta = sum(abs(delta),2);
    [J,sidx] = sort(delta);    
    match = otipData(sidx(2:11),:);    
    
    u = mean(match,1);
    s = std(match,1,1);
    e
    %{
    plot(otipData(e,:),'b')
    hold on
    errorbar(u,s,'k');
    hold off 
    drawnow
    pause(.3)
    %}
    
    e
    PRE = [PRE;u];
    plot(diag(corr(PRE,otipData(1:e,:))));
    drawnow
    
    
end


%% analysis of spec to gravi via pls
close all
X = sC;
Y = tC;
X = zscore(specData);
[Y,mu,sigma] = zscore(tipData);
[XL,YL,XS,YS,BETA,PCTVAR] = plsregress(X,Y,10);
M = [ones(size(X,1),1) X]*BETA;
M = bsxfun(@times,sigma.^-1,M);
M = bsxfun(@plus,M,mu);
%{
for e = 1:size(Y,2)
    figure;
    plot(Yfit(:,e),Y(:,e),'.');
    title([num2str(e) '--' num2str(corr(Yfit(:,e),Y(:,e)))]);
end
M = PCA_BKPROJ(Yfit,tE,tU);
delta = (M - tipData);
delta = sum(delta.*delta,2).^.5;
%}


figure;
for e = 1:size(M,1)
    plot(180/pi*M(e,:),'r');
    hold on
    plot(180/pi*tipData(e,1:DS:end),'b');
    %plot(mean(tipData,1)*180/pi,'k')
    hold off
    axis([0 61 -20 90]);
    pause(.3)
end

%{
toScore = 1;
figure;
plot(XS(:,toScore),YS(:,toScore),'.')
title(num2str(corr(XS(:,toScore),YS(:,toScore))))
%}



%% nnet
otipData = bsxfun(@minus,tipData,tipData(:,1));
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(otipData(:,1:DS:end),5);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(specData,15);
X = sC;
Y = tC;
% network
net = cascadeforwardnet([5]);
net = feedforwardnet([10 3 10]);
net.trainParam.showWindow = false;
net = fitnet(5);    
net = train(net,X',Y');
netPrediction = net(X')';
for e = 1:size(Y,2)
    figure;
    plot(netPrediction(:,e),Y(:,e),'.');
    title(num2str(corr(netPrediction(:,e),Y(:,e))));
end

M = PCA_BKPROJ(netPrediction,tE,tU);


figure;
for e = 1:size(M,1)
    plot(180/pi*M(e,:),'r');
    hold on
    plot(180/pi*(G(e,:)),'c');
    hold on
    plot(180/pi*otipData(e,1:DS:end),'b');
    plot(180/pi*tS(e,:),'k');
    %plot(mean(tipData,1)*180/pi,'k')
    hold off
    axis([0 61 -20 90]);
    pause(.3)
end

%%
figure;
for e = 1:size(M,1)
    plot(180/pi*M(e,:),'r');
    hold on
    plot(180/pi*otipData(e,1:DS:end),'b');
    %plot(mean(tipData,1)*180/pi,'k')
    hold off
    axis([0 61 -20 60]);
    pause(.3)
end


    
%% run select data sets to debug
idx = 2; % hmm fail
idx = 4; % hmm fails
idx = 1; % level set fail
idx = 100; %hmm fail long root
idx = 120; % hmm fail second root
idx = 121; % little water on tip
idx = 122; % root cap
idx = 123; % root too short
idx = 124; % root too short
idx = 125; % hmm fail on water bubble
idx = 300; % disk size
idx = 301; % size disk
idx = 302; % image flaw
idx = 304;
plot(rmData(idx,:)*180/pi);
waitforbuttonpress
load(badMat{idx},'numberFrames','numberSeedlings','CROPBOX');
SEED = badSeedlingNumber(idx);

% load the first image
[pth,nm,ext] = fileparts(badMat{idx});
imagePath = strrep(nm,'SPACE',' ');

imagePath = strrep(imagePath,'SLASH',filesep); 
imagePath = strrep(imagePath,'IA/','IA-');
imagePath = strrep(imagePath,[filesep ')'],'-)');
imagePath = [imagePath filesep];    
FileList = {};
FileExt = {'TIF','tif'};
verbose = 1;
iFileList = gdig(imagePath,FileList,FileExt,verbose);
I = imread(iFileList{1});
imshow(I,[]);
drawnow
waitforbuttonpress
curveLabelsString = ['seedling' num2str(SEED) 'curveLabels'];
fW = load(badMat{idx},curveLabelsString);
for tm = 1:numberFrames
    
    
    
    tipVecString = ['seedling' num2str(SEED) 'tipVec'];
    curveString = ['seedling' num2str(SEED) 'curve' num2str(tm)];
    tipIndexString = ['seedling' num2str(SEED) 'tipIndex'];
    f = load(badMat{idx},curveString,tipVecString,tipIndexString);
    f.(curveString).data = bsxfun(@plus,[200 50]'+CROPBOX{SEED}(1:2)',f.(curveString).data);
  
    I = imread(iFileList{tm});
    imshow(I,[]);
    hold on
    plot(f.(curveString).data(1,:),f.(curveString).data(2,:),'r');
    hold on                    
    UQg = unique(fW.(curveLabelsString){tm});
    for g = 1:numel(UQg)
        fidx = find(fW.(curveLabelsString){tm}==UQg(g));
        plot(f.(curveString).data(1,fidx),f.(curveString).data(2,fidx),[CL{UQg(g)} '*']);
        quiver(f.(curveString).data(1,f.(tipIndexString)(tm)),f.(curveString).data(2,f.(tipIndexString)(tm)),f.(tipVecString)(1,tm),f.(tipVecString)(2,tm),30)
    end
    hold off
    axis equal
    drawnow

    
end
%%
cropMaizeImage(iFileList,'',[],0,0,'/mnt/spaldingdata/nate/labelData3.mat')
%% load NAM from database
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
q = ['SELECT * FROM population_lines ' ...
                'JOIN kernel_plates ' ...
                'ON population_lines.id=kernel_plates.population_line_id ' ... 
                'JOIN kernels ' ...
                'ON kernel_plates.id=kernels.plate_id ' ...                
                'JOIN averageweightspectra_vw ' ...
                'ON kernels.id=averageweightspectra_vw.kernel_id '];                
cursor = exec(conn, q);
cursor = fetch(cursor);
fieldString = columnnames(cursor,1);
results = cursor.Data;
%% split data NAM
clear f
fidx = find(strcmp(fieldString,'plate_position'));
f.position = results(:,fidx);
fidx = find(strcmp(fieldString,'plate_name'));
f.plateName = results(:,fidx);
rm = [];
f.specData = [];
for e = 1:size(results,1)
    try
        tmp = cell2mat(results(e,35:end-1));
        f.specData = [f.specData;tmp];
    catch
        rm = [rm e];
    end
    e
end
%% random fun
[S C U E L ERR LAM] = PCA_FIT_FULL(specData,size(specData,2));
%corr(C,YS)
%%
v = double((mvnrnd(zeros(1,size(LAM,1)),diag(diag(LAM)),size(C,1))));
%v(:,256:end) = C(:,256:end);
%v = C;
synData = PCA_BKPROJ(v,E,U);
%% THIS ONE ALL WHOLE
close all
P = randperm(size(synData,1));
X = synData;
rX = specData;
nX = specData(P,:);
Y = otipData;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,10);
Ypredict = [ones(size(X,1),1),X]*BETA;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(rX,Y,10);
rYpredict = [ones(size(rX,1),1),rX]*BETA;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(nX,Y,10);
nYpredict = [ones(size(nX,1),1),nX]*BETA;
figure;plot(diag(corr(Y,Ypredict)),'b*');
hold on
plot(diag(corr(Y,rYpredict)),'r');
plot(diag(corr(Y,nYpredict)),'g');

csvwrite('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/tipAnglePredictions/factors_10_vs_rnd.csv',[diag(corr(Y,Ypredict))' diag(corr(Y,rYpredict))' diag(corr(Y,nYpredict))']);

%% show corr power ONE
close all
Y = otipData;
N = round(1000/3);
%[J sidx] = sort(rYpredict(:,50));
[J sidx] = sort(prediction(:,1));
UT = 180/pi*mean(Y(sidx(1:N),:));
ST = 180/pi*std(Y(sidx(1:N),:))/10;
UB = 180/pi*mean(Y(sidx(end-N:end),:));
SB = 180/pi*std(Y(sidx(end-N:end),:))/10;
errorbar(UT,ST,'b')
hold on
errorbar(UB,SB,'r')
%% show high oil
N = round(100/3);
[J sidx] = sort(rYpredict(:,50));
UT = mean(rYpredict(sidx(1:N),:));
ST = std(rYpredict(sidx(1:N),:))/10;
UB = mean(rYpredict(sidx(end-N:end),:));
SB = std(rYpredict(sidx(end-N:end),:))/10;
errorbar(UT,ST,'b')
hold on
errorbar(UB,SB,'r')


%% pls rand fun
rMRSE_train= [];
rMRSE_test = [];
MRSE_test_itr = [];
MRSE_train_itr = [];
rPVAL = [];
close all
X = synData;
Y = otipData;
%X = specData;

h1 = figure;
h2 = figure;
options = statset('UseParallel','always');
  

for factors = 1:50
    func = @(xtrain,ytrain,xtest,ytest)myFUN(factors,xtrain,ytrain,xtest,ytest);
    tmpTrain = {};
    tmpTest = {};
    parfor itr = 1:10
        tic
        tmp = crossval(func,X,Y,'holdout',.20);
        tmp = reshape(tmp,[size(X,1) 2]);
        fidx0 = find(tmp(:,1)==0);
        fidx1 = find(tmp(:,1)==1);
        MRSE_train_itr(itr) = mean(tmp(fidx0,2));
        MRSE_test_itr(itr) = mean(tmp(fidx1,2));
        tmpTrain{itr} = tmp(fidx0,2)';
        tmpTest{itr} = tmp(fidx1,2)';
        [h(itr) p(itr)] = ttest2(tmp(fidx0,2),tmp(fidx1,2));
        toc
    end
    rMRSE_train(:,factors) = [mean(MRSE_train_itr(:));std(MRSE_train_itr(:))];%.*numel(MRSE_train_itr).^-.5];
    rMRSE_test(:,factors) = [mean(MRSE_test_itr(:));std(MRSE_test_itr(:))];%.*numel(MRSE_test_itr).^-.5];
    for itr = 1:size(tmpTrain)
        tmpTrain{itr} = mean(tmpTrain{itr});
        tmpTest{itr} = mean(tmpTest{itr});
    end
    [h rPVAL(factors)] = ttest2(cell2mat(tmpTrain),cell2mat(tmpTest),[],[],'unequal');
    %PVAL(:,factors) = [mean(p);std(p)];
    figure(h1);
    errorbar(rMRSE_train(1,:),rMRSE_train(2,:),'r')
    hold on
    errorbar(rMRSE_test(1,:),rMRSE_test(2,:),'b')
    %errorbar(fG(1,:),fG(2,:),'g')
        
end
%% perform calibration for oil, etc
[D] = readtext('/mnt/spaldingdata/nate/121213_NIRCalibration_Spectra_NIRTube2_5reps_edited_AVG.csv');
X = cell2mat(D(2:end,26:end));
Y = D(2:end,[15:25]);
toRM = [];
for e = 1:size(Y,1)
    toRM(e) = 0;
    for j = 1:size(Y,2)
        toRM(e) = toRM(e) | isempty(Y{e,j});
    end
end
ridx = find(toRM);
Y(ridx,:) = [];
X(ridx,:) = [];
Y = cell2mat(Y);
%Y(:,7) = Y(:,1).*Y(:,9).^-1;
BETA = [];
for e = 1:size(Y,2)
    [XL,YL,XS,YS,BETA{e},PCTVAR,MSE] = plsregress(X,Y(:,e),15);
end
%predictVectors = PCA_BKPROJ(BETA(2:end,:)',caliE,caliU);
predictVectors = BETA(1:end,:)';
BETA = cell2mat(BETA);
preTest = [ones(size(prediction,1),1) specData]*BETA;
%preTest = [ones(size(X,1),1) X]*BETA;


%% find the NAM
kp = zeros(size(D,1),1);
for e =1:size(D,1)
    if ~isempty(strfind(D{e,14},'NAM'))
        if isempty(strfind(D{e,12},'Hp301'))
            if isempty(strfind(D{e,12},'P39'))
                if isempty(strfind(D{e,12},'IL14H'))
                    kp(e) =1;
                end
            end
        end
    end
    
end
fidx = find(kp);
%% test 
[ones(size(X),1) X]*BETA
%% raw rnd
close all
P = randperm(size(prediction,1));

v = double((mvnrnd(zeros(1,size(LAM,1)),diag(diag(LAM).^.5),size(C,1))));
v(:,1:100) = C(:,1:100);
synData = PCA_BKPROJ(v,E,U);

myP = [ones(size(specData,1),1) specData]*BETA;
SmyP = [ones(size(synData,1),1) synData]*BETA;


%[raw_corr raw_pre] = corr(prediction(P,:),otipData);
%[raw_corr raw_pre] = corr(prediction,otipData);
%[raw_corr raw_pre] = corr(synData,otipData);
%[raw_corr raw_pre] = corr(myP,otipData);
[raw_corr raw_pre] = corr(SmyP,otipData);
figure;
plot(raw_corr'.*(raw_pre' < .01));
figure;
plot(raw_corr')
hold on;
CL = {'r' 'b' 'g' 'k' 'm' 'c' 'r' 'b' 'g' 'k' 'm'};
figure;
for e = 1:size(raw_corr,1)
    fidx = (raw_pre(e,:) < .01);
    R = regionprops(fidx,'PixelIdxList');
    for n = 1:numel(R)
        plot(R(n).PixelIdxList,raw_corr(e,R(n).PixelIdxList),CL{e});
        hold on;
    end
end
%% corr with shovel
shD = readtext('/home/nate/Downloads/2013WiDivAngleData.csv');
shD2 = readtext('/home/nate/Downloads/WiDiv_RawRootPhenotypes_SouthAfrica_10_11_12_13.csv','[,\t]');
RS = readtext('/home/nate/Downloads/Seedling phenotyping WIDIV.csv');

%%
gpl = {};
gidx = [];
for e = 1:size(pl,2)
    tnm = pl{e};
    tnm = strrep(tnm,'-','/');
    tnm = strrep(tnm,'_','/');
    for r = 1:size(RS,1)
        if strcmp(tnm,RS{r,4})
            gpl{end+1} = RS{r,3};
            gidx = [gidx;e];
        end
    end
    e
end
%% 
MEA = [];
gMEA = [];
GNM = {};
UNM = {};
for e = 1:numel(gpl)
    tnm = gpl{e};
    tnm = upper(strrep(tnm,' ',''));
    flag = 0;
    for r = 1:size(shD,1)
        if strcmp(tnm,upper(strrep(shD{r,1},' ','')))
            tmpM = shD2(r,6:end);
            
            for l = 1:numel(tmpM)
                if strcmp(tmpM{l},'NA')
                    tmpM{l} = -1;
                end
            end
            gMEA = [gMEA;otipData(gidx(e),:)];
            MEA = [MEA;cell2mat(tmpM)];
            GNM{end+1,1} = gpl{e};
            GNM{end,2} = shD{r,1};
            UNM{end+1} = gpl{e};
            flag = 1;
            break;
        end
    end
    if flag == 0
        tnm
    end
    %e
end
ridx = find(MEA==-1);
MEA(ridx) = NaN;
%%
%BKMEA = gMEA;
%gMEA = gMEA(randperm(size(gMEA,1)),:);
UQ = unique(UNM);
UG = [];
UO = [];
for u = 1:numel(UQ)
    fidx = find(strcmp(UNM,UQ{u}));
    UG = [UG;nanmean(gMEA(fidx,:),1)];
    UO = [UO;nanmean(MEA(fidx,:),1)];
end
%UG = UG(randperm(size(UG,1)),:);
close all
LEG = {};
for e = 1:size(MEA,2)
    tmpSH = UO(:,e);
    ridx = ~isnan(tmpSH);    
    [h f] = corr(UG(ridx,:),tmpSH(ridx));
    if any(abs(h) > .1)
    
        plot(h)
        hold all
        LEG{end+1} = [shD2{1,5+e} '--' num2str(numel(find(ridx)))];
    end
    
    %LEG{1} = 'BA';
    %LEG{2} = 'CA';
end
legend(LEG);
xlabel('Time (frames)');
ylabel('Correlation with Tip Angle');
%%
close all
selidx = 1;
Y = UO(:,selidx);
kidx = ~isnan(Y);
Y = Y(kidx);
X = UG(kidx,:);
    
uTrain = [];
uTest = [];
sTrain = [];
sTest = [];
for f = 1:30    
    func = @(xtrain,ytrain,xtest,ytest)myFUN(f,xtrain,ytrain,xtest,ytest);
    MRSE_train_itr = [];
    MRSE_test_itr = [];
    for itr = 1:100    
        tmp = crossval(func,X,Y,'holdout',.20);
        tmp = reshape(tmp,[size(X,1) 2]);
        fidx0 = find(tmp(:,1)==0);
        fidx1 = find(tmp(:,1)==1);
        MRSE_train_itr(itr) = mean(tmp(fidx0,2));
        MRSE_test_itr(itr) = mean(tmp(fidx1,2));
    end
    uTrain(f) = mean(MRSE_train_itr);
    sTrain(f) = std(MRSE_train_itr);
    uTest(f) = mean(MRSE_test_itr);
    sTest(f) = std(MRSE_test_itr);
    errorbar(uTrain,sTrain,'k');
    errorbar(uTest,sTest,'r');
    hold on
    drawnow
end
%%
close all
[XL,YL,XS,YS,BETA] = plsregress(X,Y,4);
pre = [ones(size(X,1),1),X]*BETA;
plot(Y,pre,'.')
%%
close all
selidx = 2;
Y = UO(:,selidx);
kidx = ~isnan(Y);
Y = Y(kidx);
X = UG(kidx,:);

NG = 3;
kidx = kmeans(Y,NG);
%kidx = kmeans(grData(:,end),NG);
KL = [];
for u = 1:NG
    KL(u) = mean(sum(Y(kidx==u,end),2));
end
[KL sidx] = sort(KL,'descend');
kkidx = [];
for u = 1:NG
    fidx = find(kidx==u);
    kkidx(fidx) = find(sidx == u);
end
kidx = kkidx;
CMtest = [];
CMtrain = [];

CMtest = [];
[sS X sU sE sL sERR sLAM] = PCA_FIT_FULL(X,4);
%X = X(:,30:40);
for pull = 1:500
    pull
    
    trainNUM = 200;
    rndIdx = randperm(size(Y,1));
    sidx0 = rndIdx(1:trainNUM);
    sidx1 = rndIdx(trainNUM+1:end);
    
    trainX = X(sidx0,:);
    trainY = kidx(sidx0)';
    testX = X(sidx1,:);            
    testY = kidx(sidx1)';
    
    testLabel = [];
    cl = ClassificationDiscriminant.fit(trainX,trainY);    
    testLabel(:,1) = predict(cl,testX);    
    
    
    testLabel(:,2) = classify(testX,trainX,trainY);
    
    mdl = ClassificationKNN.fit(trainX,trainY,'NumNeighbors',7);
    [testLabel(:,3),score] = predict(mdl,testX);
    
    nbGau = NaiveBayes.fit(trainX, trainY);
    testLabel(:,4) = nbGau.predict(testX);
    
    tree = ClassificationTree.fit(trainX,trainY);
    testLabel(:,5) = predict(tree,testX);
    
    %testLabel = testLabel(:,3);
    testLabel = mode(testLabel,2);
    %testLabel = round(mean(testLabel,3));
    
    
    currCM = confusionmat(testY,testLabel);
    CMtest = cat(3,CMtest,currCM); %/ numel(gtestY));
    
    uJ = mean(CMtest,3)
end
%%
MAS = [];
MASS = [];
UQ = unique(pl);
for u = 1:numel(UQ)    
    Y = tipData(strcmp(UQ{u},pl),:);
    if size(Y,1) > 3
        [sS smY sU sE sL sERR sLAM] = PCA_FIT_FULL(Y,3);
        if flag = 0
            flag = 1;
            
        end
        MAS = [MAS;smY];
        MASS = [MASS;specData(strcmp(UQ{u},pl),:)];
    end
end
%% jeffs question
close all
for e = 1:100
    sidx = randperm(size(prediction,1));
    K(e) = corr(prediction(sidx(1:300),1),prediction(sidx(1:300),3));
end
hist(K)
%% theta by genotype

close all
CP = [50 250];
%CP = [64 228];
%CP = [50 200];
FX = [prediction(:,1),prediction(:,3).*prediction(:,1)];
FX = bsxfun(@minus,FX,CP);
%FX = zscore(FX);
FX(:,3) = FX(:,2).*FX(:,1).^-1;
figure;plot(FX(:,1),FX(:,2),'.');
%{
FX = [];
FX = randn([10000 2]);
FX = random('Normal',0,10,1000,2);
CP = min(FX,[],1);
FX = bsxfun(@minus,FX,CP);
FX(:,3) = FX(:,2).*FX(:,1).^-1;
%}







[~,sidx] = sort(FX(:,3));
uFX = FX;
FX = FX(sidx,:);

ABSF = cumsum(FX,1);


figure;
rmidx = isinf(FX(:,3));
FX(rmidx,:) = [];
plot(cumsum(FX(:,1)),cumsum(FX(:,2)),'.');
[p,S,u] = polyfit(cumsum(FX(:,1)),cumsum(FX(:,2)),2);
%[p] = polyfit(cumsum(FX(:,1)),cumsum(FX(:,3)),2);
%p(2) = 0;
MX = linspace(min(ABSF(:,1)),max(ABSF(:,1)),100);
MY = polyval(p,MX,S,u);
hold on
plot(MX,MY,'r');
MYi = polyval(p,ABSF(:,1),S,u);
plot(ABSF(:,1),MYi,'k.');
err = abs(MYi-ABSF(:,2));
figure;

fidx = (err < 4*10^4);
plot(ABSF(:,1),err,'m.')
hold on;
%plot(ABSF(fidx,1),err(fidx),'b.')
%plot(FX(:,1),FX(:,3),'.')
p1 = polyder(p,1);
sMY = polyval(p1,MX);
figure;
plot(MX,sMY,'k');
NP = size(FX,1);
MSX = random('Normal',mean(MX),std(MX),NP,1);
%MSXi = (MSX - u(1))/u(2);
MSSLOPE = polyval(p1,ABSF(:,1),S,u);%*mean(FX(:,2));
%MSX = random('Normal',mean(FX(:,1)),std(FX(:,1))/1.1,NP,1);

%{
MSX = random('Uniform',min(MX(1)),max(MX(end)),NP,1);
MSSLOPE = polyval(p1,MSX)*std(FX(:,2));
MSX = random('Uniform',min(FX(:,1)),max(FX(:,1)),NP,1);
%}

MSY = MSSLOPE.*FX(:,1)/u(2);%/1000;%*u(2);
%MSY = zscore(MSY);
%MSY = MSY + FX(:,2);
%MSY = MSY *std(FX(:,2));

figure;
plot(FX(:,1),FX(:,2),'.')
hold on
plot(FX(:,1),MSY,'r.');
plot(FX(fidx,1),FX(fidx,2),'k.');
%%
close all
IDX = 1:size(FX,1);
[~,IDX] = sort(sidx);
fidx = fidx(IDX);
NEWPH = ABSF(IDX,1:3);
NEWCH = FX(IDX,:);
figure;plot(corr(NEWPH,otipData)');
figure;plot(corr(NEWPH,gradient(otipData))');
newOIL = polyval(p1,NEWPH(:,1),S,u)/u(2);
figure;plot(newOIL,prediction(:,3),'.')
figure;plot(newOIL(fidx),prediction(fidx,3),'.')
figure;plot(prediction(:,1),newOIL.*prediction(:,1),'.');
figure;plot(prediction(fidx,1),newOIL(fidx).*prediction(fidx,1),'.');
figure;plot(prediction(:,1),prediction(:,1).*prediction(:,3),'.');
%%
UQ = unique(pl);
figure;
plot(NEWPH(:,1),NEWPH(:,3),'b.');
hold on
for u = 1:numel(pl)
    fidx = find(strcmp(pl,UQ{u}));
    scatter(mean(NEWPH(fidx,1)),mean(NEWPH(fidx,3)));
    hold on
    drawnow
    pause(.1)
end
%%
[x, fval, exitflag, output] = fminunc(@(p)myFITforOilWeight(p,O,W),x0);
%%
close all
CP = [50 250/100];
nPrediction = prediction;
nPrediction(:,2) = nPrediction(:,2).*nPrediction(:,1)/100;
nPrediction(:,end+1) = nPrediction(:,3).*nPrediction(:,1)/100; % 10 is mg oil
nPrediction(:,[1 10]) = bsxfun(@minus,nPrediction(:,[1 10]),CP); % subtract off min values
nPrediction(:,end+1) = nPrediction(:,10).*nPrediction(:,1).^-1*100; % 11 is percent increase oil
plot(nPrediction(:,1),nPrediction(:,10),'.')
waitforbuttonpress
toR = 11;
toRem = [toR 3];
[nPrediction,mu,sigma] = zscore(nPrediction);
XP = nPrediction(:,setdiff(1:size(nPrediction,2),toRem));
[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(XP,nPrediction(:,toR),7);
YI = [ones(size(nPrediction,1),1) XP]*BETA;
close all
plot(YI,nPrediction(:,toR),'.');
%%
waitforbuttonpress
nPrediction = bsxfun(@plus,bsxfun(@times,nPrediction,sigma),mu);
YI = YI*sigma(toR) + mu(toR);
p1 = polyfit(YI,nPrediction(:,toR),1);
close all
plot(nPrediction(:,1),nPrediction(:,10),'r.');
hold on
plot(nPrediction(:,1),YI.*nPrediction(:,1)/100,'.');
figure;
plot(YI,.5*YI.^2,'k.');
figure;
plot(corr(YI.^2,otipData),'r');
hold on
plot(corr(prediction(:,3),otipData))
NMS = {'weight','protein','oil','total_density','bone_density','bone_volume','total_volume','Starch_per','Starch_mg','mg_oil','percent_oil_increase'};
NMS = NMS(setdiff(1:size(nPrediction,2),toRem)); 
%%
subplot(size(prediction,2),size(prediction,2),1)
cnt = 1;
for e = 1:size(prediction,2)
    for e2 = 1:size(prediction,2)
        subplot(size(prediction,2),size(prediction,2),cnt)
        plot(prediction(:,e),prediction(:,e2),'.')
        cnt = cnt + 1;
    end
end
%% regress against list
[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(prediction,linspace(,8);
YI = [ones(size(prediction,1),1) XP]*BETA;

NMS = {'weight','protein','oil','total_density','bone_density','bone_volume','total_volume','Starch_per','Starch_mg'};
%% first prediction plots for shotgun create 1 for bCC
close all
plot(prediction(:,1),(prediction(:,2).*prediction(:,1)+prediction(:,3).*prediction(:,1)+prediction(:,8).*prediction(:,1))/100,'.');
hold on
plot(linspace(min(prediction(:,1)),max(prediction(:,1)),2),linspace(min(prediction(:,1)),max(prediction(:,1)),2),'r')
waitforbuttonpress
close all

cvec = ones(1,3)/3;
evec = cvec/norm(cvec);
BCC = prediction(:,[2 3 8])/100;
plot3(BCC(:,1),BCC(:,2),BCC(:,3),'.');
%{
BCC = bsxfun(@minus,BCC,cvec);
plot3(BCC(:,1),BCC(:,2),BCC(:,3),'.');
DELTA = (evec'*(BCC*evec')')';
BCC = BCC - DELTA;
BCC = bsxfun(@plus,BCC,cvec);
%}
hold on
plot3(BCC(:,1),BCC(:,2),BCC(:,3),'r.');
quiver3(cvec(1),cvec(2),cvec(3),evec(1),evec(2),evec(3),1);
waitforbuttonpress
axis([-1 1 -1 1 -1 1]);
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(BCC,3);

close all
plot(pC(:,1),pC(:,2),'.');
waitforbuttonpress
mBCC = bsxfun(@times,BCC,prediction(:,1));
mBCC(:,end+1) = prediction(:,1) - sum(mBCC,2);% 4 water fiber
mBCC(:,end+1) = prediction(:,1);                % 5 weight
mBCC(:,end+1:end+3) = bsxfun(@times,mBCC(:,1:3),(mBCC(:,end) - mBCC(:,end-1)).^-1); % percent at 6 7 8 - POS
%mBCC(:,end+1:end+3) = bsxfun(@minus,mBCC(:,1:3),(mBCC(:,end) - mBCC(:,end-1)));
plot3(pC(:,1),pC(:,2),prediction(:,1),'r.')
close all
plot3(BCC(:,1)*mean(prediction(:,1)),BCC(:,2)*mean(prediction(:,1)),BCC(:,3)*mean(prediction(:,1)),'r.')
hold on;
plot3(mBCC(:,1),mBCC(:,2),mBCC(:,3),'.');

%quiver3(cvec(1),cvec(2),cvec(3),evec(1),evec(2),evec(3),100)
quiver3(0,0,0,mean(BCC(:,1))*mean(prediction(:,1)),mean(BCC(:,2))*mean(prediction(:,1)),mean(BCC(:,3))*mean(prediction(:,1)),3);
%%
figure;
plot3(mBCC(fidx,end-2),mBCC(fidx,end-1),mBCC(fidx,end),'.');
figure;
plot3(mBCC(fidx,1),mBCC(fidx,2),mBCC(fidx,3),'.');
dL = sum(mBCC(:,1:3).*mBCC(:,1:3),2).^-.5;
mBCC(:,end+1:end+3) = bsxfun(@times,mBCC(:,1:3),dL); % 9 10 11 are normalized unit vectors
mBCC(:,end+1) = prediction(:,7); % 12 volume
mBCC(:,end+1) = prediction(:,4); % 13 total density
mBCC(:,end+1) = prediction(:,1); % 14 mass
mBCC(:,end+1) = prediction(:,1).*prediction(:,7).^-1; % 15 calc density

figure;
plot3(mBCC(:,end-2),mBCC(:,end-1),mBCC(:,end),'r.');
hold on
plot3(mBCC(:,6),mBCC(:,7),mBCC(:,8),'.');
figure;
plot3(mBCC(:,9).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,10).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,11).*(mBCC(:,5)-mBCC(:,4)),'r.');
%plot3(mBCC(fidx,9).*sum(mBCC(fidx,4:5),2),mBCC(fidx,10).*sum(mBCC(fidx,4:5),2),mBCC(fidx,11).*sum(mBCC(fidx,4:5),2),'.');
dL = sum(mBCC(:,9:11).*mBCC(:,9:11),2).^-.5;
%%
close all
figure;
plot3(mBCC(:,6),mBCC(:,7),mBCC(:,8),'b.');
hold on
plot3(mBCC(:,9).*(bindVec((mBCC(:,5)-mBCC(:,4)))+2),mBCC(:,10).*(bindVec((mBCC(:,5)-mBCC(:,4)))+2),mBCC(:,11).*(bindVec((mBCC(:,5)-mBCC(:,4)))+2),'r.')
CR = mean(mBCC(:,[6 7 8]));
quiver3(0,0,0,CR(1),CR(2),CR(3),4)


%% shotgunA
close all
figure;
SHOT = [mBCC(:,9).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,10).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,11).*(mBCC(:,5)-mBCC(:,4))];
SHOT = [mBCC(:,1),mBCC(:,2),mBCC(:,3)];
SHOT = [preTest(:,5)/100.*preTest(:,1) preTest(:,4)/100.*preTest(:,1) mBCC(:,4)];
myW = prediction(:,[1 7]); % mass volume
RATIO = mBCC(:,3).*mBCC(:,2).^-1;

rmidx = find(strcmp(genoType,'Hp301') | strcmp(genoType,'Il14H') |  strcmp(genoType,'P39') | strcmp(genoType,'bt1/+ in W23'));
rmidx = find(strcmp(genoType,'Hp301') | strcmp(genoType,'Il14H') |  strcmp(genoType,'P39') | strcmp(genoType,'bt1/+ in W23') | strcmp(genoType,' bt2/+ in W64A') | strcmp(genoType,' bt1/+ in W64A') | strcmp(genoType,' bt2/+ in W23A'));
%rmidx = prediction(:,1) <100;
%rmidx = RATIO > 700;
%rmidx = RATIO < 50;
rmidx = RATIO > 1000 | RATIO < 7
SHOT(rmidx,:) = [];
subTIP = otipData;
subTIP(rmidx,:) = [];
myW(rmidx,:) = [];

[shS shC shU shE shL shERR shLAM] = PCA_FIT_FULL(SHOT,3);
plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.')
hold on
%quiver3(shU(1),shU(2),shU(3),shE(1,1),shE(2,1),shE(3,1),200)
%quiver3(shU(1),shU(2),shU(3),shE(1,1),-shE(2,1),-shE(3,1),200)
%quiver3(shU(1),shU(2),shU(3),shU(1),shU(2),shU(3),2)
quiver3(shU(1),shU(2),shU(3),shU(1),shU(2),shU(3),2)
quiver3(shU(1),shU(2),shU(3),-shU(1),-shU(2),-shU(3),2)
shotVec = shU/norm(shU);
uSHOT = bsxfun(@minus,SHOT,shU);
deltaalong = uSHOT*shotVec';
totalong = uSHOT*shotVec' + norm(shU);
aSHOT = uSHOT - deltaalong*shotVec;
adist = sum(aSHOT.*aSHOT,2).^.5;
figure;
hist(adist,100)
DA = linspace(min(deltaalong),max(deltaalong) - .3*max(deltaalong),100);
da = 5;
h1 = figure;
h2 = figure;
h3 = figure;
plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.');
hold on
vara = [];
ua = [];
Mani = [];
for e =1:numel(DA)    
    fidx = find((deltaalong < DA(e) + da) & (deltaalong > DA(e) - da));
    subA = aSHOT(fidx,:);
    
    subMean = mean(subA,1);
    
    %subA = bsxfun(@minus,subA,subMean);
    
    subA = bsxfun(@plus,subA,shotVec*DA(e)+shU);    
    
    figure(h1);
    plot(subA(:,1),subA(:,2),'.')
    
    
    
    
    [subS subC subU subE subL subERR subLAM] = PCA_FIT_FULL(subA,3);
    theta = linspace(-pi,pi,50);
    dx = 2*(subLAM(1,1).^.5)*cos(theta);
    dy = 2*(subLAM(2,2).^.5)*sin(theta);
    subLAM = diag(diag(subLAM).^-.5);
    subLAM(3,3) = 0;
    
    plot(dx,dy,'r');
    DX = [dx' dy' zeros(size(dx'))];
    Mani(:,:,e) = PCA_BKPROJ(DX,subE,subU);
    
    UMAN(e,:) = subU;
    subE = subE*subLAM;
    EMAN(:,:,e) = subE;
    
    nadist = sum((subC.*subC),2).^.5;
    
    
    
    drawnow
    pause(.05)
    
    vara(e) = std(adist(fidx));
    ua(e) = mean(adist(fidx));
    
    vara(e) = std(nadist);
    ua(e) = mean(nadist);
    
    figure(h2);    
    plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.')
    hold on
    plot3(SHOT(fidx,1),SHOT(fidx,2),SHOT(fidx,3),'r.')
    hold off
    
    figure(h3);    
    plot3(Mani(:,1,e),Mani(:,2,e),Mani(:,3,e),'k')
    
    drawnow
end
figure(h3)
Mani = permute(Mani,[3 2 1]);
for e = 1:size(Mani,3)
    plot3(Mani(:,1,e),Mani(:,2,e),Mani(:,3,e),'k')
end
figure;
plot(vara);
figure;
plot(ua);
%%
figure;
plot3(UMAN(:,1),UMAN(:,2),UMAN(:,3),'r')
hold on
for e = 1:size(EMAN,3)
    quiver3(UMAN(e,1),UMAN(e,2),UMAN(e,3),EMAN(1,1,e),EMAN(2,1,e),EMAN(3,1,e),5)
    quiver3(UMAN(e,1),UMAN(e,2),UMAN(e,3),EMAN(1,2,e),EMAN(2,2,e),EMAN(3,2,e),5)
end
%% project data to manifold
dU = diff(UMAN,1,1);
dU = cumsum([0;sum(dU.^2,2).^.5]);
for e = 1:size(SHOT,1)
    delta = bsxfun(@minus,UMAN,SHOT(e,:));
    delta = sum(delta.*delta,2);
    [~,midx] = min(delta);
    mySHOT(e,:) = PCA_REPROJ(SHOT(e,:),EMAN(:,:,midx),UMAN(midx,:));
    mySHOT(e,3) = dU(midx);
end
figure;plot3(mySHOT(:,1),mySHOT(:,2),mySHOT(:,3),'.');
%%
[XL,YL,XS,YS,BETA] = plsregress(mySHOT,subTIP,3);
myPRE = [ones(size(mySHOT,1),1) mySHOT]*BETA;
figure;plot(myPRE')
figure;plot(corr(mySHOT,subTIP)')
%%


close all
toP = [1 2 3];
plot3(mBCC(:,toP(1)),mBCC(:,toP(2)),mBCC(:,toP(3)),'.');
figure;
plot(prediction(:,1),mBCC(:,4),'.')
plot(prediction(:,1),mBCC(:,1),'.')
plot(prediction(:,1),mBCC(:,2),'.')
plot(prediction(:,1),mBCC(:,3),'.')
plot3(prediction(:,1),mBCC(:,toP(1)),mBCC(:,toP(3)),'.');
figure;
plot(mBCC(:,1),mBCC(:,4),'.')

%%
[x1,x2] = ndgrid(linspace(1,10,100),linspace(1,10,100));
f = x2.*x1.^-1;
mesh(x1,x2,f);
%%
close all
[x1,x2] = ndgrid(linspace(-1,1,50),linspace(-1,1,50));
f = x2.*x1.^-1;
mesh(x1,x2,f);
%%
close all
toM = 3;
toD = 1;
oil = prediction(:,toD).*prediction(:,toM)/100;
[x1,x2] = ndgrid(linspace(min(prediction(:,toD)),max(prediction(:,toD)),20),linspace(min(oil),max(oil),20));
f = x2.*x1.^-1;
gf1 = -x2.*x1.^-2;
gf2 = x1.^-1;

df1 = @(x1,x2)-x2.*x1.^-2;
df2 = @(x1,x2)x1.^-1;
dl = @(x1,x2)(df1(x1,x2)*df1(x1,x2) + df2(x1,x2)*df2(x1,x2)).^-.5;
df1n = @(x1,x2)df1(x1,x2)*dl(x1,x2);
df2n = @(x1,x2)df2(x1,x2)*dl(x1,x2);
mesh(x1,x2,f);

%[gf1 gf2] = gradient(f);
dl1 = (gf1.*gf1 + gf2.*gf2).^-.5;
gf1 = gf1.*dl1;
gf2 = gf2.*dl1;
dl1 = (gf1.*gf1 + gf2.*gf2).^-.5;
initP = [];
for l = 1:size(x1,1)    
    %initP = [x1(1,1) x2(1,1)];
    %initP = [175 150];
    initP(1,:,l) = [x1(l,1) x2(l,1)];
    initP(1,:,l) = [x1(1,l) x2(1,l)];
    dt = 1;
    for e = 1:100
        dx1 = df1n(initP(e,1,l),initP(e,2,l))*dt;
        dx2 = df2n(initP(e,1,l),initP(e,2,l))*dt;
        dx1 = -df2n(initP(e,1),initP(e,2))*dt;
        dx2 = df1n(initP(e,1),initP(e,2))*dt;
        %initP = [initP;initP(e,:) + [dx1 dx2]];
        initP(e+1,:,l) = initP(e,:,l) - [dx1 dx2];
    end
end


hold on
plot3(prediction(:,toD),oil,oil.*prediction(:,toD).^-1,'r.');
%quiver(x1,x2,gf1,gf2);
%quiver(x1,x2,-gf1,gf2);
for l = 1:size(initP,3)
    plot(initP(:,1,l),initP(:,2,l),'k')
end
axis([min(x1(:)) max(x1(:)) min(x2(:)) max(x2(:)) min(oil)/100 max(oil)/100])
figure;
contour(x1,x2,f,100)
hold on
plot(prediction(:,toD),oil,'r.','MarkerSize',1)
plot(initP(:,1),initP(:,2),'k')
%quiver(x1,x2,gf2,gf1);
%quiver(x1,x2,-gf1,gf2);
%%
close all
toM = 3;
toD = 1;
oil = prediction(:,toD).*prediction(:,toM)/100;
starch = prediction(:,toD).*prediction(:,8)/100;
[x1,x2] = ndgrid(linspace(min(oil),max(oil),20),linspace(min(starch),max(starch),20));
f = x1 + x2;
mesh(x1,x2,f);
hold on
plot3(oil,starch,oil+starch,'r.');
%%
axis([min(x1(:)) max(x1(:)) min(x2(:)) max(x2(:)) min(oil)/100 max(oil)/100])
figure;
%%
close all
contour(x1,x2,f,50)
hold on
plot(oil,starch,'k.','MarkerSize',2)
%%
%% shotgunB
close all
figure;
SHOT = [mBCC(:,9).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,10).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,11).*(mBCC(:,5)-mBCC(:,4))];
SHOT = [mBCC(:,1),mBCC(:,2),mBCC(:,3),mBCC(:,4)];
myW = prediction(:,[1 7]);
RATIO = mBCC(:,3).*mBCC(:,2).^-1;

rmidx = find(strcmp(genoType,'Hp301') | strcmp(genoType,'Il14H') |  strcmp(genoType,'P39') | strcmp(genoType,'bt1/+ in W23'));
rmidx = find(strcmp(genoType,'Hp301') | strcmp(genoType,'Il14H') |  strcmp(genoType,'P39') | strcmp(genoType,'bt1/+ in W23') | strcmp(genoType,' bt2/+ in W64A') | strcmp(genoType,' bt1/+ in W64A') | strcmp(genoType,' bt2/+ in W23A'));
%rmidx = prediction(:,1) <100;
%rmidx = RATIO > 700;
%rmidx = RATIO < 50;
rmidx = RATIO > 1000 | RATIO < 7
SHOT(rmidx,:) = [];
subTIP = otipData;
subTIP(rmidx,:) = [];
myW(rmidx,:) = [];

[shS shC shU shE shL shERR shLAM] = PCA_FIT_FULL(SHOT,4);
shotVec = shU/norm(shU);
uSHOT = bsxfun(@minus,SHOT,shU);
deltaalong = uSHOT*shotVec';
totalong = uSHOT*shotVec' + norm(shU);
aSHOT = uSHOT - deltaalong*shotVec;
adist = sum(aSHOT.*aSHOT,2).^.5;
figure;
hist(adist,100)
DA = linspace(min(deltaalong),max(deltaalong) - .3*max(deltaalong),100);
da = 25;
h1 = figure;
h2 = figure;
h3 = figure;
plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.');
hold on
vara = [];
ua = [];
Mani = [];
for e =1:numel(DA)    
    fidx = find((deltaalong < DA(e) + da) & (deltaalong > DA(e) - da));
    subA = aSHOT(fidx,:);
    
    subMean = mean(subA,1);
    
    %subA = bsxfun(@minus,subA,subMean);
    
    subA = bsxfun(@plus,subA,shotVec*DA(e)+shU);    
    
    figure(h1); 
    
    
    
    
    [subS subC subU subE subL subERR subLAM] = PCA_FIT_FULL(subA,4);
    theta = linspace(-pi,pi,50);
    dx = 2*(subLAM(1,1).^.5)*cos(theta);
    dy = 2*(subLAM(2,2).^.5)*sin(theta);
    subLAM = diag(diag(subLAM).^-.5);
    subLAM(3,3) = 0;
    
    plot(dx,dy,'r');
    DX = [dx' dy' zeros(size(dx'))];
    Mani(:,:,e) = PCA_BKPROJ(DX,subE,subU);
    
    UMAN(e,:) = subU;
    subE = subE*subLAM;
    EMAN(:,:,e) = subE;
    
    nadist = sum((subC.*subC),2).^.5;
    
    
    
    drawnow
    pause(.05)
    
    vara(e) = std(adist(fidx));
    ua(e) = mean(adist(fidx));
    
    vara(e) = std(nadist);
    ua(e) = mean(nadist);
    
    figure(h2);    
    plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.')
    hold on
    plot3(SHOT(fidx,1),SHOT(fidx,2),SHOT(fidx,3),'r.')
    hold off
    
    figure(h3);    
    plot3(Mani(:,1,e),Mani(:,2,e),Mani(:,3,e),'k')
    
    drawnow
end
figure(h3)
Mani = permute(Mani,[3 2 1]);
for e = 1:size(Mani,3)
    plot3(Mani(:,1,e),Mani(:,2,e),Mani(:,3,e),'k')
end
figure;
plot(vara);
figure;
plot(ua);


%%
HELLOW = [];
for e = 1:size(pl,2)
    if ~isempty(strfind(pl{e},'07S'))
        HELLOW = [HELLOW;[prediction(e,[1 7 4]) preTest(e,[1 9 7])]];
        isNAM(e) = 1;
    else
        isNAM(e) = 0;
    end
    e
end
HELLOW = [];
for e = 1:size(pl,2)
    if ~isempty(strfind(pl{e},'WISN'))
        HELLOW = [HELLOW;[prediction(e,[1 7 4]) preTest(e,[1 9 7])]];
        isWISN(e) = 1;
    else
        isWISN(e) = 0;
    end
    e
end
figure;
plot(HELLOW(:,1).*HELLOW(:,2).^-1,HELLOW(:,3),'.')
hold on
plot(HELLOW(:,1).*HELLOW(:,2).^-1,HELLOW(:,6),'r.')
plot(HELLOW(:,4).*HELLOW(:,5).^-1,HELLOW(:,6),'k.')
%%
plot(HELLOW(:,2),HELLOW(:,4),'k.')
%%
close all
[nS nC nU nE nL nERR nLAM] = PCA_FIT_FULL(specData,5);
TOTS = sum(specData,2);
DF = linspace(min(TOTS),max(TOTS),50);
dS = std(TOTS)/2;
newU = PCA_REPROJ(2*nU,nE,nU);
AO = PCA_REPROJ(ones(size(nU)),nE,nU);
%%

for e = 1:numel(DF)
    fidx = find((TOTS < DF(e) + dS) & (TOTS > DF(e) - dS));
    plot3(nC(:,1),nC(:,2),nC(:,3),'b.');
    hold on
    plot3(nC(fidx,1),nC(fidx,2),nC(fidx,3),'r.')
    quiver3(0,0,0,newU(1),newU(2),newU(3))
    quiver3(0,0,0,AO(1),AO(2),AO(3),'g')
    drawnow
    pause(.2)
end
%{
%% test with gravitropsim
[nS nC nU nE nL nERR nLAM] = PCA_FIT_FULL(otipData,3);
TOTS = sum(otipData,2);
NV = PCA_REPROJ(otipData,nU'/norm(nU),nU);
figure;
plot3(NV,TOTS,1*ones(size(NV,1),1),'b.')
hold on
plot3(NV(logical(isNAM)),TOTS(logical(isNAM)),2*ones(sum(isNAM),1),'r.')
plot3(NV(logical(isWISN)),TOTS(logical(isWISN)),3*ones(sum(isWISN),1),'k.')
plot3(NV(logical(~isWISN & ~isNAM)),TOTS(logical(~isWISN & ~isNAM)),4*ones(sum(~isWISN & ~isNAM),1),'c.')
%}
myVAL = corr(specData,nU'/norm(nU));
%%
NV = PCA_REPROJ(specData,nU'/norm(nU),nU);
figure;
plot3(NV,TOTS,1*ones(size(NV,1),1),'b.')
hold on
plot3(NV(logical(isNAM)),TOTS(logical(isNAM)),2*ones(sum(isNAM),1),'r.')
plot3(NV(logical(isWISN)),TOTS(logical(isWISN)),3*ones(sum(isWISN),1),'k.')
plot3(NV(logical(~isWISN & ~isNAM)),TOTS(logical(~isWISN & ~isNAM)),4*ones(sum(~isWISN & ~isNAM),1),'c.')
%%
close all
specX = specData(logical(isNAM),:);
[nS nC nU nE nL nERR nLAM] = PCA_FIT_FULL(specX,5);
TOTS = sum(specX,2);
DF = linspace(min(TOTS),max(TOTS),50);
dS = std(TOTS)/2;
newU = PCA_REPROJ(2*nU,nE,nU);
%%
for e = 1:numel(DF)
    fidx = find((TOTS < DF(e) + dS) & (TOTS > DF(e) - dS));
    plot3(nC(:,1),nC(:,2),nC(:,3),'b.');
    hold on
    plot3(nC(fidx,1),nC(fidx,2),nC(fidx,3),'r.')
    quiver3(0,0,0,newU(1),newU(2),newU(3))
    drawnow
    pause(.2)
end
%%
NV = PCA_REPROJ(specX,nU'/norm(nU),nU);
figure;
plot3(NV,TOTS,ones(size(NV,1),1),'b.')
%% try total engergy in
close all
mySPEC = bsxfun(@times,specData(find(isWISN),:),sum(specData(find(isWISN),:),2).^-1);
[nS nC nU nE nL nERR nLAM] = PCA_FIT_FULL(mySPEC,5);
plot3(nC(:,1),nC(:,2),nC(:,3),'b.');
%%
close all
X = specData(:,450:end);
X = specData(:,1:450);
X = specData;
%X = diff(specData,2,2);
Y = otipData;
%[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,200);
[S Y U Ex L ERR LAM] = PCA_FIT_FULL(Y,3);
%X = X - 1;
%X = exp(X).^-1;
%X = bsxfun(@times,X,mean(X,1) < 1) + bsxfun(@times,-(X-1) + 1,mean(X,1) > 1);

%X = zscore(X);
%Y = zscore(Y);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    %set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end

for e = 1:size(mV,2)
    figure;
    bar(mB(:,e));
    %set(gca,'XTickLabel',LAB.ICA,'XTick',1:numel(LAB.ICA));
    title(['Coff Shape:' num2str(e)]);
end


[Ycore YcoreP] = corr(mU,otipData);
%[Ycore YcoreP] = corr(mU,Y);
for e = 1:size(mB,2)
    figure;
    %bar(Ycore(e,:))
    plot((Ycore(e,:)));
    %set(gca,'XTickLabel',LAB.shape(1:3),'XTick',1:numel(LAB.shape(1:3)))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
%[Xcore XcoreP] = corr(mU,preTest(pidx,subIDX));
[Xcore XcoreP] = corr(mU,preTest);
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    %plot(bar(Xcore(e,:)));
    set(gca,'XTickLabel',NMS,'XTick',1:numel(NMS))
    %set(gca,'XTickLabel',NMS(subIDX),'XTick',1:numel(subIDX))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(corr(mU(:,e),mV(:,e)))]);
end

%[S Y U Ex L ERR LAM] = PCA_FIT_FULL(otipData,3);
[XL,YL,XS,YS,beta,PCTVAR] = plsregress(mU,Y,3);
preTip = [ones(size(mU,1),1) mU]*beta;
preTip = PCA_BKPROJ(preTip,Ex,U);
figure;plot(preTip');
figure;plot(diag(corr(preTip,otipData)));

%plot3(mU(:,1),mU(:,2),mU(:,3),'.')
%csvwrite('~/T.csv',[mean(specData,1)' mA]);
%% investigate
close all
plot(mean(X,1))
dY = zeros(size(specData,1),15,size(specData,2));
for e =1:size(specData,1)
    dY(e,:,:) = cwt(specData(e,:),1:15,'gaus2');
    e
end
%%
close all
xi = linspace(0,1,1000);
yi = log(xi);
plot(xi,yi)
%%
hold on
plot(bindVec(mA(:,1)) + min(mean(specData,1)) - mean(bindVec(mA(:,1))),'r')
figure;
plot(bindVec(mA(:,1)))
hold on
plot(mean(specData,1) < 1,'r')
%%
X = specData;
Y = otipData;

%[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,30);
[S Y U Ex L ERR LAM] = PCA_FIT_FULL(Y,3);
%X = zscore(X);
%Y = zscore(Y);
[XL,YL,XS,YS,beta,PCTVAR] = plsregress(X,Y,3);
preTip = [ones(size(X,1),1) X]*beta;
preTip = PCA_BKPROJ(preTip,Ex,U);
figure;plot(preTip');
figure;plot(diag(corr(preTip,otipData)));
%% hold out for something
close all

%X = preTest2(:,subIDX);
%X = preTest;
%Y = [T.shapeDataL(:,1:3)];

%Y = grData(:,end);


%[S Y U Ex L ERR LAM] = PCA_FIT_FULL(Y,3);


%Y = bsxfun(@minus,Y,mean(Y,1));
%X = bsxfun(@minus,X,mean(X,1));
%{
for e = 1:size(Y,1)
    Y(e,:) = Y(e,:) * norm(Y(e,:))^-1;
end
for e = 1:size(X,1)
    X(e,:) = X(e,:) * norm(X(e,:))^-1;
end
%}
UQ = unique(pl);

fG = [];
fG2 = [];
fG3 = [];

MASFG = [];
MASFG2 = [];
MAD = [];
UMT = [];
UMB = [];
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
h5 = figure;
MC = [];
for factors = 3:600
    X = specData;    
    [S X U Ex L ERR LAM] = PCA_FIT_FULL(X,factors);
    Y = otipData;
    [S Y U Ex L ERR LAM] = PCA_FIT_FULL(Y,3);
    for pull = 1:50
        fG = [];
        fG2 = [];
        fG3 = [];
        R = randperm(numel(UQ));
        R = R(1:50);
        sidx0 = [];
        sidx1 = [];
        for e = 1:numel(R)
            %sidx0 = [sidx0 find(~strcmp(pl,UQ{R(e)}))];
            sidx1 = [sidx1 find(strcmp(pl,UQ{R(e)}))];
        end
        sidx0 = setdiff(1:numel(pl),sidx1);
        trainX = X(sidx0,:);
        trainY = Y(sidx0,:);
        testX = X(sidx1,:);        
        %testX = mean(testX,1);
        testY = Y(sidx1,:);
        %testY = mean(testY,1);

        uSx = mean(trainX);
        uSy = mean(trainY);
    
        [A,B,r,trainX,trainY,stats] = canoncorr(trainX,trainY);
        testX = bsxfun(@minus,testX,uSx);
        testY = bsxfun(@minus,testY,uSy);
        testX = testX*A;
        testY = testY*B;
        [tC] = corr(testX,testY);
        if any(isnan(tC(:)))
            report = 1;
            holVx = testX;
            holVy = testY;
            
            break;
        end
        MC(pull,:,factors-2) = diag(tC);
        MCMOD(pull,:,factors-2) = r';
    end
    factors;
    uMC = mean(MC,1);
    uMCMOD = mean(MCMOD,1);
    plot(squeeze(uMC)')
    hold on
    plot(squeeze(uMCMOD)','--')
    hold off
    drawnow
end