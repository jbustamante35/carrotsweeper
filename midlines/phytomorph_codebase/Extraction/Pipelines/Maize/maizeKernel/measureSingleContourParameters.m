function [M] = measureSingleContourParameters(contour,B)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                measureSingleContourParameters.m  (Inputs are relative to 1200dpi)
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
        currCentroid = fliplr(mean(contour));
        VEC = [fliplr(contour(1,:)) - (currCentroid)];
        VEC = 3*VEC;
        H = norm(VEC);

        % create line to sample along both major and minor                        
        nVEC = [-VEC(2) VEC(1)];

        % make line to sample along
        nX = linspace(currCentroid(1)-nVEC(1),currCentroid(1)+nVEC(1),2*H);
        nY = linspace(currCentroid(2)-nVEC(2),currCentroid(2)+nVEC(2),2*H);

        % sample along the binary image and measure the major axis
        nLP = ba_interp2(double(B),nX,nY);
        nLPf = imfill(~logical(nLP>.5),[1 round(numel(nLP)/2)]);
        nfidx = find(nLP.*nLPf>.5);

        % make line to sample along
        X = linspace(currCentroid(1)-VEC(1),currCentroid(1)+VEC(1),2*H);
        Y = linspace(currCentroid(2)-VEC(2),currCentroid(2)+VEC(2),2*H);                

        % sample along the binary image and measure minor axis
        LP = ba_interp2(double(B),X,Y);
        LPf = imfill(~logical(LP>.5),[1 round(numel(LP)/2)]);
        fidx = find(LP.*LPf>.5);

        % plot the major and minor axis
        M.MajorLine = [[nX(nfidx(1));nX(nfidx(end))],[nY(nfidx(1));nY(nfidx(end))]];
        M.MinorLine = [[X(fidx(1));X(fidx(end))],[Y(fidx(1));Y(fidx(end))]];        

        % stack the major and minor axis
        M.MajorLength = numel(nfidx)-1;
        M.MinorLength = numel(fidx)-1;

        % get the eth contour and project to the kernels frame
        % of reference
        tmp = interp1(linspace(0,1,size(contour,1)),contour,linspace(0,1,500));
        utmp = mean(tmp,1);
        VEC = VEC/norm(VEC);
        nVEC = nVEC/norm(nVEC);
        C = PCA_REPROJ(tmp,[VEC',nVEC']',utmp);
        M.kernelReferenceFrame = [VEC',nVEC']';
        M.kernelMean = utmp;
        M.iContour = C; 
        M.contour = contour;
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:measureSingleContourParameters.m******\n']);
    end
end