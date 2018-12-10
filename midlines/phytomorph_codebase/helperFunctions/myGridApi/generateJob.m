function [jP] = generateJob(para,Ostore)
    try
        %import phytoG.locked.Bcompute.resultsContainer;        
        %resultsContainer = resultsContainer();
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % note: here generate job can be a few ideas
        % the job may have been externally generated.
        % in this case generateJob should load the jobInformation
        % and return result. in other case, generate job should generate job     
        % and return result.
        %%%%%%%%%%%%%%%%%%%%%%%%    

        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % configure@inPort@common
        %%%%%%%%%%%%%%%%%%%%%%%%
        inPara.language    = 'matlab';                
        inPara.fileList    = {};
        inPara.fileExt     = {'tiff','TIF','PNG','png','tif'};
        inPara.returnType  = 'set';
        inPara.verbose     = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % configure@outPort@common
        %%%%%%%%%%%%%%%%%%%%%%%%        
        para.localFileBase = [];
        para.request_uniqueKey = 1;
        para.basePath = [para.basePath filesep para.jobType filesep datestr(now) filesep];
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % configure@job set@common
        %%%%%%%%%%%%%%%%%%%%%%%%
        jP = myHS_X('pJob');
        switch para.jobType
            
            case 'Seed Size'
                % assume that para.inPort.fileList is a cell array of files
                %%%%%%%%%%%%%%%%%%%%%%%%
                % configure@inPort@common and indiv
                %%%%%%%%%%%%%%%%%%%%%%%%
                iport = inPort();
                iport.sourceType = 'fileList';
                iport.source = para.fileList;
                % pull from inPort
                iport.pull();
                iport.fileList.xForm_iRods(para.iPlantUser);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % configure@outPort@indiv                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                oport = outPort();
                oport.setWriteBase(para.basePath);                
                oport.setUser(para.iPlantUser);                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % configure@jobObject@indiv
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                jP{1} = pJob();
                jP{1}.setIport(iport);
                jP{1}.setOport(oport);
                jP{1}.attachFunction(@(x)seedSize_main(x));
                jP{1}.setDisp(1);
            case 'Root Morphometrics'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate job for each in the fileList
                for e = 1:numel(para.fileList)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@inPort@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    inPara.filePath = {[para.fileList{e} filesep]};
                    iport = inPort();
                    iport.sourceType = 'structuredPathScan';
                    iport.source = inPara;                    
                    % pull from inPort                    
                    iport.pull();
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@outPort@indiv                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [tp,tnm,text] = fileparts(iport.fileList{1}.fileName);
                    oport = outPort();
                    oport.setWriteBase(para.basePath);
                    oport.genUniqueKey(tp);
                    oport.setUser(para.iPlantUser);                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@jobObject@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    jP{e} = pJob();
                    jP{e}.setIport(iport);
                    jP{e}.setOport(oport);
                    jP{e}.attachFunction(@(x)stackProcess(x));
                    jP{e}.setDisp(1);
                end
                % gather draw data
                jP.executeStagingFunction();
            case 'Root Morphometrics V2'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate job for each in the fileList
                for e = 1:numel(para.fileList)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@inPort@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    inPara.filePath = {[para.fileList{e} filesep]};
                    iport = inPort();
                    iport.sourceType = 'structuredPathScan';
                    iport.source = inPara;                    
                    % pull from inPort                    
                    iport.pull();
                    % transform filenames to irods file names
                    iport.fileList.xForm_iRods(para.iPlantUser);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@outPort@indiv                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [tp,tnm,text] = fileparts(iport.fileList{1}.fileName);
                    oport = outPort();
                    oport.setWriteBase(para.basePath);
                    oport.genUniqueKey(tp);
                    oport.setUser(para.iPlantUser);                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@jobObject@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    jP{e} = pJob();
                    jP{e}.setIport(iport);
                    jP{e}.setOport(oport);
                    jP{e}.attachFunction(@(x)processImageStack(x));
                    jP{e}.setDisp(1);
                    %jP{e}.setResultsContainer(resultsContainer);
                end
                % gather draw data
                %jP.executeStagingFunction();
            case 'Kinematics with Steady State Flow'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate job for each in the fileList
                for e = 1:numel(para.fileList)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@inPort@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    inPara.filePath = {[para.fileList{e} filesep]};
                    iport = inPort();
                    iport.sourceType = 'structuredPathScan';
                    iport.source = inPara;                    
                    % pull names from inPort
                    iport.pull();
                    % transform filenames to irods file names
                    iport.fileList.xForm_iRods(para.iPlantUser);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@outPort@indiv                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [tp,tnm,text] = fileparts(iport.fileList{1}.fileName);
                    oport = outPort();
                    oport.setWriteBase(para.basePath);
                    oport.genUniqueKey(tp);
                    oport.setUser(para.iPlantUser);                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@jobObject@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    jP{e} = pJob();
                    jP{e}.setIport(iport);
                    jP{e}.setOport(oport);
                    jP{e}.attachFunction(@(x)straightKinematics_main(x));
                    jP{e}.attachStagingFunction(@(x)gather_drawingData(x));
                    jP{e}.attachParaField('kinematicsType ','steadyState');
                    jP{e}.setDisp(1);                    
                end
                % gather draw data
                jP.executeStagingFunction();
                % issue itickets over the image files
                %jP.issueTickets(2);
                %jP.postJob(para.iPlantUser);
                %jP.postJob('phytotest');
            case 'Kinematics without Steady State Flow'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate job for each in the fileList
                for e = 1:numel(para.fileList)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@inPort@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    inPara.filePath = {[para.fileList{e} filesep]};
                    iport = inPort();
                    iport.sourceType = 'structuredPathScan';
                    iport.source = inPara;                    
                    % pull from inPort                    
                    iport.pull();
                    iport.fileList.xForm_iRods(para.iPlantUser);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@outPort@indiv                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [tp,tnm,text] = fileparts(iport.fileList{1}.fileName);
                    oport = outPort();
                    oport.setWriteBase(para.basePath);
                    oport.genUniqueKey(tp);
                    oport.setUser(para.iPlantUser);                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@jobObject@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    jP{e} = pJob();
                    jP{e}.setIport(iport);
                    jP{e}.setOport(oport);
                    jP{e}.attachFunction(@(x)straightKinematics_main(x));
                    jP{e}.attachStagingFunction(@(x)gather_drawingData(x));
                    jP{e}.attachParaField('kinematicsType ','nonSteadyState');
                    jP{e}.setDisp(1);
                end
                % gather draw data
                jP.executeStagingFunction();
                % issue itickets over the image files
                jP.issueTickets(2);
                %jP.postJob(para.iPlantUser);
                jP.postJob('phytotest');
            case 'Kinematics newI'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate job for each in the fileList
                for e = 1:numel(para.fileList)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@inPort@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    inPara.filePath = {[para.fileList{e} filesep]};
                    iport = inPort();
                    iport.sourceType = 'structuredPathScan';
                    iport.source = inPara;                    
                    % pull from inPort                    
                    iport.pull();
                    %if isdeployed
                        iport.fileList.xForm_iRods(para.iPlantUser);
                    %end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@outPort@indiv                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [tp,tnm,text] = fileparts(iport.fileList{1}.fileName);
                    oport = outPort();
                    oport.setWriteBase(para.basePath);
                    oport.genUniqueKey(tp);
                    oport.setUser(para.iPlantUser);                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@jobObject@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    jP{e} = pJob();
                    jP{e}.setIport(iport);
                    jP{e}.setOport(oport);
                    jP{e}.attachFunction(@(x)straightKinematics_main_newI(x));
                    jP{e}.attachStagingFunction(@(x)gather_drawingData(x));
                    jP{e}.attachParaField('kinematicsType ','nonSteadyState_nonLogistic');
                    jP{e}.setDisp(1);
                end
                % gather draw data
                jP.executeStagingFunction();
                %if isdeployed
                    % issue itickets over the image files
                    jP.issueTickets(2);
                    %jP.postJob(para.iPlantUser);
                    %jP.postJob('phytotest');
                %end
                
            case 'Tip Tracker'
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate job for each in the fileList
                for e = 1:numel(para.fileList)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@inPort@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    inPara.filePath = {[para.fileList{e} filesep]};
                    iport = inPort();
                    iport.sourceType = 'structuredPathScan';
                    iport.source = inPara;                    
                    % pull from inPort                    
                    iport.pull();
                    iport.fileList.xForm_iRods(para.iPlantUser);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@outPort@indiv                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [tp,tnm,text] = fileparts(iport.fileList{1}.fileName);
                    oport = outPort();
                    oport.setWriteBase(para.basePath);
                    oport.genUniqueKey(tp);
                    oport.setUser(para.iPlantUser);                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % configure@jobObject@indiv
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                    jP{e} = pJob();
                    jP{e}.setIport(iport);
                    jP{e}.setOport(oport);
                    jP{e}.attachFunction(@(x)tip_tracking_main(x));
                    jP{e}.attachStagingFunction(@(x)gather_drawingData_tipTracking(x));
                    jP{e}.setDisp(1);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % gather draw data
                jP.executeStagingFunction();
            case 'Train Particle Learner'
                %%%%%%%%%%%%%%%%%%%%%%%%
                % configure@outPort@special
                %%%%%%%%%%%%%%%%%%%%%%%% 
                rmidx = strfind(para.basePath,filesep);
                para.basePath(rmidx(end-1)+1:rmidx(end)) = [];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % configure@inPort@indiv
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                inPara.filePath = {[para.basePath filesep para.jobType]};
                iport = inPort();
                iport.sourceType = 'structuredPathScan';
                iport.source = inPara;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % configure@outPort@indiv                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                oport = outPort();
                oport.setWriteBase(para.basePath);                
                oport.setUser(para.iPlantUser); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % configure@jobObject@indiv
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sessionFile = [para.basePath 'session.mat'];
                if ~exist(sessionFile)
                    jP{1} = pJob();
                    jP{1}.setIport(iport);
                    jP{1}.setOport(oport);
                    jP{1}.attachFunction(@(x)learnGUI(x));
                    %jP{1}.attachFunction(@(x)generate_learning_set(x));
                    %jP{1}.attachStagingFunction(@(x)gather_learnerData(x));
                    jP{1}.setDisp(1);
                    jP{1}.attachParaField('localMountPoint',para.local_irodsMNT);
                    jP{1}.attachParaField('iPlantUser',para.iPlantUser);
                    jP{1}.attachParaField('basePath',para.basePath);
                else
                    data = load(sessionFile);
                    jP = data.d;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % gather draw data
                jP.executeStagingFunction();
            case 'mongoDB'
                import phytoG.locked.Bpersist.Bos.implementations.*;
                iport = inPort();
                iport.sourceType = 'mongoDB';
                iport.source = OStore_mdb();
                iport.id = '4f617fa83d1e9d1666715faa';
                iport.pull();
        end    
    catch ME
        getReport(ME)
    end
end