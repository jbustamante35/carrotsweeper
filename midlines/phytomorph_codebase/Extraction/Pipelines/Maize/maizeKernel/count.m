function [fidx sel uF] = count(feature)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                count.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                feature:      The information is needed. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        sr = [];
        for k = 1:numel(feature)
            RAT = round(feature.*feature(k).^-1);
            sr(k) = sum(RAT==1);    
        end
        MO = mode(sr);
        fidx = (sr == MO);
        uF = mean(feature(fidx));
        sel = round(feature.*uF^-1);    
        fidx = (sel==1);
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:count.m******\n']);
    end
end