function [jobList,fileList] = convert2BSONPtr(jobList,fileList)
    import phytoCompute.*;
    % each set
    for s = 1:numel(fileList)
        
        % load job or make new jobPackage
        if jobList.containsField(num2str(s))            
            jbPackage = jobList.get(jb);
        else
            jbPackage = jobPackage();
        end
        
        % each file : add to var list and return ptr to value
        for e = 1:numel(fileList{s})
            fileList{s}{e} = jbPackage.addFileVar(fileList{s}{e});
        end
        
        % add job package
        jobList.put(num2str(s),jbPackage);
    end
end 