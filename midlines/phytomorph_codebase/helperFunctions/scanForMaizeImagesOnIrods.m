function [FileList] = scanForMaizeImagesOnIrods(user,plantType,tissueType)
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
    end
    FileList = {};
    FileExt = {'tiff','TIF','tif'};
    verbose = 1;
    FileList = gdig(pathToScan,FileList,FileExt,verbose);
end