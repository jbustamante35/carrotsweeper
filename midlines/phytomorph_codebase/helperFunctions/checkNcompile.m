function [] = checkNcompile(mexaName,curPath)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                checkNcompile.m checks existance of .mexa64 file and if it
                does not then compile and create .mexa64.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                cwtK_closed_imfilter.m, circshift.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                dB:      The information is needed. 
                B:
                I:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        % hard coded destination
        mexaName = ['/mnt/snapper/Lee/gitHub_maizepipeline/maizePipeline/helperFunctions/ba_interp/' mexaName];
        if exist(mexaName, 'file')
            % File exists.  Do stuff....
            cd(curPath);
        else
          % File does not exist.
          %cd('/mnt/snapper/Lee/gitHub_maizepipeline/maizePipeline/helperFunctions/ba_interp');
          mex -O /mnt/snapper/Lee/gitHub_maizepipeline/maizePipeline/helperFunctions/ba_interp/ba_interp2.cpp;
          %mex -O ba_interp2.cpp;
          addpath(mexaName);
          cd(curPath);
        end
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:checkNcompile.m******\n']); 
    end
end