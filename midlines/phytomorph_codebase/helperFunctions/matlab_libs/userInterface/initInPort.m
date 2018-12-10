function [job] = initInPort(varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select if not selected
    if nargin == 0
        selectedPaths = uipickfiles();
    else 
        selectedPaths = varargin{1};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select outport
    outDir = uigetdir();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct job
    for e = 1:numel(selectedPaths)
        % inport
        job(e).inPort.language     = 'matlab';       
        job(e).inPort.filePath     = {selectedPaths{e}};
        job(e).inPort.fileList     = {};
        job(e).inPort.fileExt      = {'tiff','TIF','png'};    
        job(e).inPort.returnType   = 'set';
        job(e).inPort.verbose      = 1;
        job(e).inPort.fileList     = structuredPathScan(port(e).inPort);
        % outport
        job(e).outPort.path        = outDir;
        job(e).outPort.specialCharater = '--';
        job(e).outPort.nDepth = numel(para.inPort.fileList{1});
    end
    
end