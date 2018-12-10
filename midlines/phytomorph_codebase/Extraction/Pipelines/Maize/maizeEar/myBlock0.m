function [ret] = myBlock0(block)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                myBlock0.m takes block from siteProcess and return mean of fft.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                bsxfun.m, fft.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                block:       An image to be analyze in a string that includes path and file name.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        % subtract off the mean
        uBlock = mean(block,1);
        block = bsxfun(@minus,block,uBlock);
        % perform fft along 1 dim
        fT = fft(block,[],1);
        % get the mean of the fft signal along the 2nd dim
        ret = mean(abs(fT),2);
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:myBlock0.m******\n']);
    end
end