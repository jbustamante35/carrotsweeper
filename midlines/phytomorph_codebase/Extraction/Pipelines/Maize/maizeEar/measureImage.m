function [T ret] = measureImage(fsI,toMeasure,downsample,dR,CHUNK)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                measureImage.m is main function to handle cob analysis. It takes all input variables 
                for its dependent functions. This function returns final result including image with 
                bounding box and color circle. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                siteProcess.m, myBlock0.m, findT.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                fsI:            The image, period to be measured.
                toMeasure:      Errode such that the fft window samples only ear image.
                downsample:     The number of down sample grid sites.
                dR:             Set the current window size
                CHUNK:          The number of chunk for input for FFT in myBlock0.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    %%%%%%%%%%%%%%%%%%%%%%%
    % init return vars    
    T = NaN;   
    ret = NaN;
    %%%%%%%%%%%%%%%%%%%%%%% 
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set display to false
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create sample sites for fft application
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [g1 g2] = ndgrid(1:downsample:size(toMeasure,1),1:downsample:size(toMeasure,2));
        idx = sub2ind(size(toMeasure),g1(:),g2(:));
        idx = find(toMeasure(idx) == 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform fft block process at sitesret{k} = func(B);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;        
        sites = [g1(idx) g2(idx)];
        ret = siteProcess(sites,[dR 1],@(block)myBlock0(block),fsI,CHUNK);
        ret = cell2mat(ret);
        ret = mean(ret,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % analysis of the average fft over the windows
        % look at the first half of the period fft signal 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ufT = ret(1:dR);
        % filter the first half of the signal
        h = fspecial('average',[5 1]);
        ufT = imfilter(ufT,h);
        % find period via first peak finding
        [T f] = findT(ufT,2*dR+1);
        fprintf(['e-time for find fft @ ' num2str(etime(clock,tm)) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % display
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if disp
            imshow(sI,[]);
            hold on
            CS = 1:size(fsI,1);
            plot(100*cos(CS*T^-1*2*pi) + 100,CS,'r');
            title(num2str(T));
            hold off
            drawnow
        end
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:measureImage.m******\n']);
    end
end