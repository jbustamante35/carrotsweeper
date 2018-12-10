function [germNet frameNet frameChain] = calibrate(calFile,matFileList)
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) load the calibration file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    calFile = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/Gustin_20170425_emergence_hand_score.csv';
    J = readtext(calFile);
    dataPaths = unique(J(2:end,1));
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.5) spin up cFlow for condor run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    func = cFlow('processEmergenceStack_v3');
    func.setMCRversion('v840');
    func.setMemory(8000);
    miniStack = {};
    miniMask = {};
    
    numToRun = numel(dataPaths);
    %numToRun = 2;
    
    for j = 1:numToRun
        fprintf(['starting spool job number:' num2str(j) '\n']);tic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2) scan for images and sort based on str2num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dataPath = dataPaths{j};
        CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
        [o,r] = system(CMD);
        [r] = parseRecords(r);  
        for e = 1:numel(r)
            FileList{e} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
            [~,n] = fileparts(FileList{e});
            nm(e) = str2num(n);
        end
        [~,sidx] = sort(nm);
        FileList = FileList(sidx);
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3) issue tickets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [FileList] = issueBulkTicket(FileList);
        % local run
        %[miniStack{j},miniMask{j}] = processEmergenceStack_v3(FileList,[]);
        % condor run
        [miniStack{j},miniMask{j}] = func(FileList,[]);
        
        fprintf(['ending spool job number:' num2str(j) ':@' num2str(toc) '\n']);
    end
    
    
    
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    func.submitDag(auth,50,50);
    
end