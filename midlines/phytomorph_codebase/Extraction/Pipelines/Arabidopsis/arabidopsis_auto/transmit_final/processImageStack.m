function [out] = processImageStack(jobStruct)
        try
        stack = jobStruct.inPort.fileList;
        %saveStruct,disp

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init display
        if jobStruct.disp
            disp = figure;
        else
            disp = [];
        end

        %%%%%%%%%%%%%%%%%%%
        % extract contours and midlines
        for e = 1:numel(stack)
            out{e} = processImage(stack{e},disp);
            fprintf(['done with image@' num2str(e) ':' num2str(numel(stack)) '\n'])
        end

        segLen = 20;
        % for each time frame
        for e = 1:numel(stack)
            % for each midline
            for m = 1:numel(out{e}.midlines)
                midline = out{e}.midlines(m).data;            
                % compute the angle
                dM = diff(midline,1,2);
                seg = mean(-dM(:,1:segLen),2);
                angle(e,m) = atan2(seg(2),seg(1))*180/pi;
                % compute the length
                dL = sum(dM.*dM,1).^.5;
                dL = sum(dL);
                length(e,m) = dL;
            end        
        end
        growthRate = gradient(length);
        % 
        jobStruct.outPort.push('midline_contour.mat',out);
        jobStruct.outPort.push('angle.csv',angle);
        jobStruct.outPort.push('length.csv',length);
        jobStruct.outPort.push('growthRate.csv',growthRate);
        %%%%%%%%%%%%%%%%%%%
        %
        %[pth,nm,ext] = fileparts(stack{1});
        %fidx = strfind(pth,filesep);
        %pth = pth(fidx(end-saveStruct.DEPTH)+1:end);
        %pth = strrep(pth,filesep,'----');

        %%%%%%%%%%%%%%%%%%%
        %
        %matOutPath = [saveStruct.outPath 'mat/'];
        %mkdir(matOutPath);
        %fileName = [matOutPath pth '.mat'];


        %%%%%%%%%%%%%%%%%%%
        %
        %spoolToDisk(out,fileName);
        catch ME
            rethrow(ME);
        end
end