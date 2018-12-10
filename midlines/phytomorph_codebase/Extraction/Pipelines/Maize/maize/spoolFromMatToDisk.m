    function [] = spoolFromMatToDisk()
    %% load from condor mini - post hmm - startLOAD here
    close all
    inFilePath = ['/mnt/snapper/nate/myDev/maizeWhole_mini6/output/'];
    inFilePath = '/mnt/spaldingdata/nate/mJUNK/';
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
    toKeep = [];
    fprintf(['Start\n']);
    %%%%%%%%%%%%%%%%%%%%%
    % for each file
    %%%%%%%%%%%%%%%%%%%%%
    for m = 1:numel(mFileList)
        m
        numel(mFileList)
        load(mFileList{m},'numberFrames','numberSeedlings','CROPBOX');
        tipAngle = [];
        stackFlag = 1;
        tmp_toKeep = [];
        try
            
            %%%%%%%%%%%%%%%%%%%%%
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
                    % load the curve for obtaining the growthrate
                        tipPath = [];
                        for tm = 1:frDS:numberFrames
                            curveString = ['seedling' num2str(s) 'midline' num2str(tm)];
                            tmpData = load(mFileList{m},curveString);
                            tmpData = tmpData.(curveString);
                            tmpData = diff(tmpData,1,2);
                            tmpData = sum(tmpData.*tmpData,1).^.5;
                            tipPath = [tipPath sum(tmpData)];
                        end
                        dP = tipPath;

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
                                    stackFlag = 0;
                                end
                            end
                        end

                        % stack the tip angle 
                        tipAngle = [tipAngle;fW.(tipAngleString)];
                        masterName = {masterName{:} names{s}};
                        setupName = {setupName{:} SETUP};
                        plateName = {plateName{:} tmpPlateName};
                        seedlingNumber = [seedlingNumber s];
                        matFile{end+1} = mFileList{m};
                        iPath{end+1} = imagePath;
                        tmp_toKeep = [tmp_toKeep;1];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % end dispi
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % catch error and display loading error
                    catch ME
                        stackFlag = 0;
                        mFileList{s};
                        disp(ME.getReport);
                    end
                    % try to stack the growth data
                    if stackFlag
                        try
                            MasterGrowthSur = [MasterGrowthSur;dP];
                        catch ME
                            stackFlag = 0;
                        end
                    end
                    if ~stackFlag
                        tmp_toKeep(:) = 0;
                        break
                    end
                end
                if stackFlag
                    plot(180/pi*tipAngle');
                    axis([0 61 -30 90]);
                    title(num2str(numberSeedlings));
                    drawnow
                    %waitforbuttonpress;
                    MasterTipAngle = [MasterTipAngle;tipAngle];
                end
                toKeep = [toKeep;tmp_toKeep];
                fprintf(['Size of mat file:' num2str(numel(matFile)) ':' num2str(size(MasterTipAngle,1)) ':' num2str(numel((toKeep))) '\n']);
            end
        catch ME
            ME
        end
    end
 %{
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
%}
%% backup mat file
matFileBK = matFile;
%% spool out for Takeshi
matFile = matFile(find(toKeep));
%%
oPath = '/mnt/piranhapuffer2/Bessie/MaizeData/';
oPath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/fullExtract_12.03.15/';
mkdir(oPath);
UQ = unique(matFile);
for u = 1:numel(UQ)
    fidx = find(strcmp(matFile,UQ{u}));
    [pth,nm,ext] = fileparts(UQ{u});
    fileName = [oPath nm '.csv'];
    fileName = strrep(fileName,'SLASH','_');
    fileName = strrep(fileName,'SPACE',' ');
    %tmpTipdata = fMA(fidx,:);
    %tmpTipdata = fMG(fidx,:);
    tmpTipdata = MasterTipAngle(fidx,:);
    csvwrite(fileName,tmpTipdata');
    %plot(tmpTipdata');
    %drawnow
    u
    numel(UQ)
end

end