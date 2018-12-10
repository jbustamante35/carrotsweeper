function [finalScore,rawX,nX] = gradeRedCaps(PCA_scores,emergenceNet,BESTEmergence,popThreshold,whenThresh)
    if ischar(PCA_scores)
        load(PCA_scores,'rawX');
        PCA_scores = rawX;
    end
    nX = {};
    rawX = PCA_scores;
    pulseResponse = [];
    finalScore = [];
    for e = 1:(size(PCA_scores,3))
        
        tmpL = PCA_scores(:,:,e);
        %tmpL = PCA_scores(:,(1:170),e);
        tsig = convertToStrain(tmpL,11,10);
        dsig = convertToDoesEmergeSignal(tmpL,11,10);


        doesEmergeProb{e} = emergenceNet.predict(dsig);
        doesEmerge{e} = emergenceNet.classify(dsig);
        doesEmerge{e} = double(double(doesEmerge{e})==1);
       
        whenEmerge{e} = BESTEmergence.predict(tsig);


        %MXemPROB = max(doesEmergeProb{e}(1),whenEmerge{e}(2,end));
        MXemPROB = doesEmergeProb{e}(1);
        doesEmerge{e} = double(MXemPROB > popThreshold);


       
        BESTEmergence.resetState();
        emergenceNet.resetState();


        [stepSig,finalScore(e)] = extractFrameFromProb(whenEmerge{e}(2,:),whenThresh,5,21);
        finalScore(e) = finalScore(e) * doesEmerge{e};
        pulseResponse(:,e) = stepSig'* doesEmerge{e};


        %{
        init = tmpL(:,1);
        tmpL = bsxfun(@minus,tmpL,init);
        tmpL = bsxfun(@times,tmpL,init.^-1);
        tmpL = imfilter(tmpL,fspecial('average',[1 11]),'replicate');
        % stack both QR and zscore QR
        %}
        nX{e} = tsig;
    end
    
end
%{
matFile = '/home/nate/Downloads/20180130_Rack2_Camera1_dataPackage.mat'
finalScore = gradeRedCaps(matFile,emergeNET,trainedNet_FM,.6);


FilePath = '/mnt/tetra/nate/RED/';
redMATFileList = {};
FileExt = {'mat'};
redMATFileList = gdig(FilePath,redMATFileList,FileExt,1);

redTHRESH = linspace(.51,.8,20);
JfinalScore = {};
oPath = '/mnt/tetra/nate/forJEFFredCapps/'
mkdir(oPath)
for t = 1:numel(redTHRESH)
    tPath = [oPath 'threshold' num2str(t) filesep];
    mkdir(tPath)
    parfor e = 1:numel(redMATFileList)
        tic
       
        [pth,nm,ext] = fileparts(redMATFileList{e});
        tName = [tPath nm '.csv'];
        [JfinalScore{t,e},~,~] = gradeRedCaps(redMATFileList{e},emergeNET,trainedNet_FM,redTHRESH(t));
        csvwrite(tName,JfinalScore{t,e});
        toc
    end
end

%}