function [resultsPath] = ipWrap(inputFile)
    try        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % read the input
        % IN.F := filelist in ID numbers
        % IN.P := filelist in ID numbers
        %myf = xml_load(inputFile);
        %%%%%%%%%%%%%%%%%%%%%%%%
        % generate matlab commands
        %cmd = gen_myf(myf);
        %%%%%%%%%%%%%%%%%%%%%%%%
        % call commands
        for i = 1:numel(myf)        
            eval(cmd{i});
        end
    catch ME
        ME
    end

end

