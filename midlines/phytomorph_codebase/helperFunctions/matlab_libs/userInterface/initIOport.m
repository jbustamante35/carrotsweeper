function [job] = initIOport(varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select if not selected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 0
        selectedPaths = uipickfiles();
    else 
        selectedPaths = varargin{1};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select outport
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outDir = uigetdir();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct job
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(selectedPaths)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inport
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        job(e).inPort.language     = 'matlab';       
        job(e).inPort.filePath     = {selectedPaths(e).pathData};
        job(e).inPort.fileList     = {};
        job(e).inPort.fileExt      = {'tiff','TIF','png','tif'};    
        job(e).inPort.returnType   = 'list';
        job(e).inPort.verbose      = 1;
        fileList = structuredPathScan(job(e).inPort);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over found sets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for s = 1:numel(fileList)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % outport
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        job(e).outPort.path        = outDir;
        job(e).outPort.specialCharater = '--';
        job(e).outPort.nDepth      = numel(para.inPort.fileList{1});
    end
    
end