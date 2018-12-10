function [objStore] = ios(iV)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % import java libs        
    import phytoG.locked.Bpersist.Bos.implementations.*;
    %%%%%%%%%%%%%%%%%%%%%%%%
    % create output ports    
    switch iV.storeFlag
        case 'dbSQL'
            mkdir(iV.outPath);
            dbFile = [iV.outPath getDateVec() '_out.db'];
            objStore = OlStore_sql(dbFile);
        case 'mongoDBlocal'
            objStore = OStore_mdb();
        case 'mongoDBremote'
            objStore = OStore_mdbR();
    end
    objStore.accessResource();
end