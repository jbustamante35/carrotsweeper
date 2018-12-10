function [FileList] = scanForImagesOnIrods(user,plantType,tissueType,FileExt)
    
    FileList = {};
    
    switch tissueType
        case 'ears'
            tissueType = 'ear';                        
        case 'cobs'
            tissueType = 'cob';
        case 'kernels'
            tissueType = 'kernel';
        case 'seedlings'
            tissueType = 'seedling';
        case 'wholes'
            tissueType = 'whole';
        case 'rootStraight'
            pathToScan = [pathToScan 'kinematics/straight/'];
            %pathToScan = [pathToScan 'kinematics/singleMidline/'];
            %pathToScan = [pathToScan 'kinematics/singleMidline/Aditi/Hypo and Hyper/'];
        case 'stomataTopo'
            tissueType = 'stomataTopo';
        case 'kernelSwellData'
            tissueType = 'kernel_swellData';
        case 'crossSection'
            tissueType = 'crosssectionData';
    end
    [r] = getRawDataList(user,plantType,tissueType);
    for e = 1:numel(r)
        [p,nm,ext] = fileparts(r(e).DATA_NAME);
        if any(strcmp(ext(2:end),FileExt))
            FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        end
    end
    
    
    %{
    CMD = gMeta(FileList{1},'add','key','value');
    CMD = uncLog(FileList(1:100),'add','kernelAnalysis','1.0',[],[],0);
    CMD = uncLog(FileList(1:100),'rmw','kernelAnalysis','1.0',[],[],0);
    CMD = uncLog(FileList(1),'add','kernelAnalysis','1.0',{'1'},{'0'},0);
    %}
    
    
    %{
    pathToScan = '/home/nate/iplant/';
    CMD = ['mountIrods.sh /iplant/home/' user '/$plantTypeData/'];
    CMD = strrep(CMD,'$plantType',plantType);
    [o,r] = system(CMD);
    switch tissueType
        case 'ears'
            pathToScan = [pathToScan 'earData/'];                        
        case 'cobs'
            pathToScan = [pathToScan 'cobData/'];
        case 'kernels'
            pathToScan = [pathToScan 'kernelData/'];
        case 'seedlings'
            pathToScan = [pathToScan 'seedlingData/'];
        case 'wholes'
            pathToScan = [pathToScan 'wholeData/'];
        case 'rootStraight'
            pathToScan = [pathToScan 'kinematics/straight/'];
            %pathToScan = [pathToScan 'kinematics/singleMidline/'];
            %pathToScan = [pathToScan 'kinematics/singleMidline/Aditi/Hypo and Hyper/'];
        case 'stomataTopo'
            pathToScan = [pathToScan 'stomataTopo/'];
        case 'kernelSwellData'
            pathToScan = [pathToScan 'kernel_swellData/'];
        case 'crossSection'
            pathToScan = [pathToScan 'crosssectionData/'];
    end
    FileList = {};
    verbose = 1;
    FileList = gdig(pathToScan,FileList,FileExt,verbose);
    %}
end