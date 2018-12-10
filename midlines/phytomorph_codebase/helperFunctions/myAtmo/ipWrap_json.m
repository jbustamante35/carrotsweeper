function [resultsPath] = ipWrap(f_handle,inputFile)
    % this function will wrap matlab scripts
    % f_handle :=  function to call
    % inputFile : = xml file for running code
    
    %%%
    % STEP 1: if the inputFile is not given
    % then assume that it is on the working path
    if nargin == 1;        
        
        cdir = dir('./');
        cdir([cdir.isdir]) = [];
        L = struct2list(cdir,'name');
        
        % find the first xml file
        for i = 1:numel(L)
            [pth nm ext] = fileparts(L{i});
            if strcmp(ext,'.xml')
                kp(i) = 1;
            else
                kp(i) = 0;
            end
        end
        inputFile = L(fidx(1));
        
    end
    
    %%%
    % read the input
    % IN.F := filelist in ID numbers
    % IN.P := filelist in ID numbers
    myf = xml_load(inputFile);
    
    %%%
    % init workspace
    workingPath = ['.' filesep 'worker_working' filesep];    
    resultsPath = myf(1).arg{2};
    %pwd
    mkdir(workingPath);
    mkdir(resultsPath);
    
    %%%
    % pre-clean local
    %fprintf(['workingPath:' workingPath '\n']);
    %pause(5);
    %localClean(workingPath);
    
    %%%
    % bring in data to working directory
    % note: if from iRODS then use icommands else pass through from local system
    myf(1).arg{1} = xfer_get(myf(1).arg{1},workingPath,0);
    
    %%%
    % call main function
    cmd = gen_myf(myf);
    for i = 1:numel(myf)
        %cmd = gen_myf(myf(i));
        eval(cmd{i});
    end
    %{
    %%%
    % put data
    xfer_put(resultsPath,remotePath)
    %}
    
    % post-clean local
    localClean(workingPath);
    
    

end

