function [] = generateImageClassifier(fileList)

    stackMasterList = {};
    imageType = {};
    
    for e = 1:numel(fileList)
        [pth,nm,ext] = fileparts(fileList{e});
        pth = [filesep pth];
        fidx = strfind(pth,filesep);
        if ~isempty(strfind(pth,'imageStacks'))
            stackMasterList{end+1} = fileList{e};
            imageType{end+1} = pth((fidx(2)+1):(fidx(3)-1));
        end
        
    end
    
    
    
    UQ = unique(imageType);
    for e = 1:numel(imageType)
        imageTypeLabel(e) = find(strcmp(imageType{e},UQ));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the image-list data source 
    imageSource = dataSource(stackMasterList,@(X)double(imread(X)),1);
    % create nozzle for dataSource
    thumb_Nozzle = dataNozzle(@(X,e0,e1)func_thumbNail(X,[50 50],1,false),imageSource,1);
    % get the reduction nozzle for the thumb_nozzle
    rThumb_nozzle = thumb_Nozzle.getReductionNozzle(3);
    % transform nozzle into transduction nozzle
    % self-note: comment the transduction nozzle
    tThumb_nozzle = rThumb_nozzle.getTransductionNozzle(@(X)getTrainedNB(X,imageTypeLabel));
    
end

%{
    FilePath = '/mnt/snapper/nate/phunny/';
    FileList = {};
    FileExt = {'tif','TIF'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    generateImageClassifier(FileList)
%}