function [ret] = structuredPathScan(inPort)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % inPort
    %               mainBase    :> location to scan
    %               language    :> matlab or java
    %               type        :> set or list
    %%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % scan inPort
    %%%%%%%%%%%%%%%%%%%%%%%%    
    switch inPort.language
        case 'matlab'
            ret = [];
            %%%%%%%%
            % scan from matlab            
            verbose = 1;
            %%%%%%%%
            % switch return type
            switch inPort.returnType
                case 'set'                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % perform set scan in matlab
                    % create image stack response variable
                    % pour data into the imageStacks
                    ret = myHS_X('imageStack');
                    for p = 1:numel(inPort.filePath)
                        % perform dig
                        tmp = sdig(inPort.filePath{p},inPort.fileList,inPort.fileExt,verbose);                                                
                        for s = 1:numel(tmp)
                            % create image stack
                            iS = imageStack();
                            % pour images into stack
                            for i = 1:numel(tmp{s})
                                iS{i} = imageFile(tmp{s}{i});
                            end
                            % sort stack
                            iS.sort();
                            % return stack
                            ret{numel(ret)+1} = iS;                            
                        end
                        
                        
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if number con then filter
                    if isfield(inPort,'numCon')
                        if ~isempty(inPort.numCon)
                            % look for sets with given number of images
                            for e = 1:numel(ret)
                                k(e) = numel(ret{e}) == inPort.numCon;
                            end
                            % return on
                            ret = ret(find(k));
                        end
                    end
                    
                    
                    
                    
                case 'list'
                    for p = 1:numel(inPort.filePath)
                        tmp = gdig(inPort.filePath{p},inPort.fileList,inPort.fileExt,verbose);
                        for r = 1:numel(tmp)
                            ret{end+1} = tmp{r};
                        end
                    end                
            end
        case 'java'
            %%%%%%%%
            % scan from java
            import phytoG.locked.Bpersist.Bfs.*;
            import phytoG.locked.Bpersist.Bfs.implementations.*;
            % init the scanner
            localFS = Bfs_local();
            scanner = fsScanner();
            % add the scan path
            scanner.addScanPath(cellTojarray(inPort.filePath,'java.lang.String'));
            scanner.addExt(cellTojarray(inPort.fileExt,'java.lang.String'));
            % init the scanner
            scanner.initialize();
            scanner.setFileSystem(localFS);
            %%%%%%%%
            % switch return type
            switch inPort.returnType
                case 'set'              
                    scanner.setScanType(javaMethod('valueOf', 'phytoG.locked.Bpersist.Bfs.fsScanner$scanType', 'foTds'));
                    if isfield(inPort,'numCon')
                        if isempty(inPort.numCon)
                            setImageNumberConstraint(inPort.numCon);
                        end
                    end
                    % set the fileSystem            
                    ret = scanner.setScan(); 
                    
                case 'list'
                    scanner.setScanType(javaMethod('valueOf', 'phytoG.locked.Bpersist.Bfs.fsScanner$scanType', 'fiTds'));
                    % set the fileSystem            
                    ret = scanner.listScan();
            end
    end
end
%{
    %%%%%%%%%%%%%%%%%%%%%%%%
    % inPort
    %               mainBase    :> location to scan
    %               language    :> matlab or java
    %               type        :> set or list
    %%%%%%%%%%%%%%%%%%%%%%%%

    inPort.language     = 'java';
    inPort.filePath     = {'/mnt/spaldingdata/genLab/algorithms/Rafael/raw_data/toProcess/'};
    inPort.fileList     = {};
    inPort.fileExt      = {'tiff','TIF'};    
    inPort.returnType   = 'list';
    inPort.verbose      = 1;

    structuredPathScan(inPort);
    %%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%
    % setup paths
    %%%%%%%%%%%%%%%%%%%%%%%%
    switch inPort.moniker
        case 'trieupham10'
            inPort.mainBase = '/mnt/scratch3/users/trieupham10/phytoMorph/morphoMetrics/';
        case 'bessie'            
            inPort.mainBase = '/mnt/spaldingdata/genLab/algorithms/Bessie/phytoMorph/morphoMetrics/';
        case 'moonsoo'
            inPort.mainBase = '/mnt/spaldingdata/Moonsoo/phytoMorph/morphoMetrics/';
        case 'rafael'
            inPort.mainBase = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/';
        case 'TaN'
            inPort.mainBase = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/';
        case 'Rue'
            inPort.mainBase = '/mnt/spaldingimages/nate/whole_rue_set/';
        case 'myTest'
            inPort.mainBase = '/mnt/spaldingdata/genLab/algorithms/Rafael/';
        case 'custom'
            inPort.mainBase = inPort.mainBase;
    end


%}
    