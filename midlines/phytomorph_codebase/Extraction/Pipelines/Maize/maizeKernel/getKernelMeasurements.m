function [M] = getKernelMeasurements(B,contours) 
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                getKernelMeasurements.m  (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                cwtK_closed_imfilter.m, circshift.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                B:      The information is needed. 
                contours:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        for k = 1:numel(contours)
            M(k) = measureSingleContourParameters(contours{k},B);
        end    
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:getKernelMeasurements.m******\n']);
    end
end