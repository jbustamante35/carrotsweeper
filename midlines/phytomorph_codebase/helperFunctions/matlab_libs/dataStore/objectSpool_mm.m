function [] = objectSpool_mm(outPort,data)
    try
        report = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % libs
        import phytoG.locked.Bpersist.Bos.implementations.*;
        import phytoG.locked.BdataObjects.BbioObjects.arabidopsis.*;        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make object store
        localOS = OStore_fs(outPort.objPath);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over each image
        for img = 1:numel(data.toDisk(data.toStore.IDX).d)
            for obj = 1:numel(data.toDisk(data.toStore.IDX).d{img})
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % create root apex model
                rApex = mmRootApex(localOS);
                for phen = 1:numel(data.toDisk(data.toStore.IDX).d{img}{obj})
                    if report
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % generate reporting string
                        rep = ['img:' num2str(img)      ...
                               '__@obj:' num2str(obj)     ...
                               '__@phen:' num2str(phen)   ...
                               '__@type:' data.toDisk(data.toStore.IDX).d{img}{obj}{phen}.class ...
                               '__@name:' data.toDisk(data.toStore.IDX).d{img}{obj}{phen}.label ...
                               '\n'];                    
                        fprintf(rep);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % generate object
                    out = mat2pct(data.toDisk(data.toStore.IDX).d{img}{obj}{phen}.data, ...
                            data.toDisk(data.toStore.IDX).d{img}{obj}{phen}.class);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % project object to disk
                    localOS.put(out);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % put into rApex model
                    if ~isempty(data.toDisk(data.toStore.IDX).d{img}{obj}{phen}.icmd)
                        fnc = ['rApex.' data.toDisk(data.toStore.IDX).d{img}{obj}{phen}.icmd '(out);'];
                        eval(fnc);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % save data model
                localOS.put(rApex);
            end
        end
    catch ME
        
    end
end