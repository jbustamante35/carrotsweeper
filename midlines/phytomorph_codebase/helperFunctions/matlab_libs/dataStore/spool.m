function [fileList] = spool(outPort,data)
    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        import phytoG.locked.Bpersist.Bos.abstraction.*;
        fileList = {};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make outport directory
        if ~isIRODS(outPort.basePath)
            mkdir(outPort.basePath);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make unique output directory
        if isfield(outPort,'uniqueKey')
            kPath = [outPort.basePath outPort.uniqueKey];            
        else
            kPath = [outPort.basePath];
        end
        if ~isIRODS(kPath);mkdir(kPath);end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write out csv files to disk
        if outPort.disk
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            kpPath = [kPath filesep 'csv_image_Data' filesep];
            if ~isIRODS(kpPath);mkdir(kpPath);end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call to write data
            fileList = myWriter(kpPath,data.toDisk);
        end
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write out objects to store
        if outPort.oStore        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make outport directory
            mkdir(outPort.basePath);                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make object directory
            objPath = [kPath filesep 'objStore' filesep];
            mkdir(objPath);            
            outPort.objPath = objPath;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call to spool
            outPort.toStore(outPort,data);
            target_ext = char(OStore.ext);
            fileList = gdig(outPort.objPath,fileList,{target_ext(2:end)},1);
        end
        %}
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % package to zip
        if outPort.zip
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % project to zip
            zipFile = [kPath filesep 'return.zip'];
            projectTozip(fileList,zipFile);
        end
        %}

    
    catch ME
        ME
    end
end