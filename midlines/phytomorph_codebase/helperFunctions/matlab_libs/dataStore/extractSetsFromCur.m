function [fileListSet firstImgSet] = extractSetsFromCur(cur,dataStore,N,randFlag)
    import phytoGraph_lite.Bpersist.Bos.implementations.OStore_mdb.*;
    import phytoConnectivity.PortForwarding.*;
    import phytoGraph_lite.BdataStructures.implementations.*;
    import phytoGraph_lite.BBioData.Bimg_dataSet.*;   
    import ch.ethz.ssh2.*;
    %%%%%%%%%%%%%%%%%%%%%%%
    % init object
    firstImgSet = java.util.ArrayList;
    fileListSet = {};
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % loop over cursor
    for imgSet = 1:N
        %%%%%%%%%%%%%%%%%%%%%%%
        % extract first dataset from cusor
        fileList = java.util.ArrayList;        
        iS = Bimg_dataSet(cur.next(),dataStore);
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % sort iS on file name
        %%%%%%%%%%%%%%%%%%%%%%%
        for img = 1:iS.size()    
            tmp = char(iS.getFile(img-1).getFileName());
            [p,n,e] = fileparts(tmp);
            fn(img) = str2num(n);
        end
        [fn sidx] = sort(fn);

        %%%%%%%%%%%%%%%%%%%%%%%
        % make the list of first images
        firstImgSet.add(iS.get(sidx(1)-1));

        %%%%%%%%%%%%%%%%%%%%%%%    
        % add to list in sorted order
        %%%%%%%%%%%%%%%%%%%%%%%
        for img = 1:iS.size()  
            fileList.add(iS.getFile(sidx(img)-1));            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % add to the first image set
        fileListSet{imgSet} = fileList;
        
    end
    
end