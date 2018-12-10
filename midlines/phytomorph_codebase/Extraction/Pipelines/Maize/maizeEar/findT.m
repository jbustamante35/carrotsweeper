function [T f Tchoices SNIP] = findT(sig,N,varargin)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                findT.m is main function to handle ear analysis. It takes all input variables 
                for its dependent functions. This function returns final result including image with 
                bounding box. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                bindVec.m, interp1.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                sig:       	
                N:           
                varargin:          Input arguments passed.
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        MAXT = 600; % max kernel length
        if nargin == 3
            MAXT = varargin{1};
        end
        Tline = N./((1:numel(sig))-1)';

        T_thresh = Tline < MAXT;
        % changed on aug 11 2016 from 8 to 10
        [localMAX] = nonmaxsuppts(sig,20);

        % added for simulation
        %sig = sig.*T_thresh;

        nsig = bindVec(sig);
        thresh = graythresh(nsig);
        bidx = (nsig > thresh);


        fidx = find(localMAX.*bidx.*T_thresh);
        fpeak = sig(fidx);
        [fpeak sidx] = sort(fpeak,'descend');

        fidx

        f = mean(fidx(sidx(1)));
        % first peak
        f = fidx(1);
        %%
        T = N/(f-1);
        f = T^-1;
        if nargout >= 3
            Tchoices = N.*fidx.^-1;
        end
        if nargout == 4
            SNIP = interp1(Tline.^-1,sig,linspace(0,.025,1000));
        end
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:findT.m******\n']);
    end
end