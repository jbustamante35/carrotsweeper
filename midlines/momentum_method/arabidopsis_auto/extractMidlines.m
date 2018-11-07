function [db] = extractMidlines(contour,tP)   
    try        
        report.time = 1;
        report.verbose = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract midlines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INPUT: 
        %           I       := image
        %           para    := parameters for running the script         
        %                   := para.sig         -> sigma for gaussian filter
        %                   := para.gradPara    -> gradient configure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OUTPUT: 
        %           cM      := corner map       -> corner strength
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load java libraries
        import phytoG.locked.BdataObjects.BbioObjects.arabidopsis.*;
        import phytoG.locked.Bpersist.Bos.implementations.*;        
        if report.time;tm = clock;end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% init for objective function for midline
        GEN = 20;
        POP = 10;
        MAG = 2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % isolate POINTS - tip and base(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%
            % label the contour tensor
            %%%%%%%%%%%%%%%%%
            H = rootLabelProcessChain();
            P = H(contour);
            P{1} = tP{1};
            %%%%%%%%%%%%%%%%%
            % upperBase point
            %%%%%%%%%%%%%%%%%
            base = P{2}{1};
            %%%%%%%%%%%%%%%%%
            % tip
            %%%%%%%%%%%%%%%%%
            tip = P{1}{1};
            %%%%%%%%%%%%%%%%%
            % lowerBase
            %%%%%%%%%%%%%%%%%
            lowerBase = P{3}{1};
            %%%%%%%%%%%%%%%%%
            % base point
            %%%%%%%%%%%%%%%%%
            upperBase = P{4}{1};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % isloate CURVES - midlines and contours
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            %%%%%%%%%%%%%%%%%
            sec_d = contour(P{4}{2}:P{3}{2},:);
            %%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%
            % splice out the upper contour segment
            %%%%%%%%%%%%%%%%%
            mag_upper = flipud(sec_d(1:P{1}{2},:));
            mag_upper = arcLength(mag_upper,'mag',MAG);
            %%%%%%%%%%%%%%%%%                
            % splice out the lower contour segment
            %%%%%%%%%%%%%%%%%
            mag_lower = (sec_d(P{1}{2}:end,:));
            mag_lower = arcLength(mag_lower,'mag',MAG);
            %%%%%%%%%%%%%%%%%
            % match boundary
            %%%%%%%%%%%%%%%%%
            matched_lower = boundaryManifold(mag_lower,mag_upper,GEN,POP,0);
            matched_upper = boundaryManifold(mag_upper,mag_lower,GEN,POP,0);
            %%%%%%%%%%%%%%%%%        
            % get midlines
            %%%%%%%%%%%%%%%%%
            midline1 = (matched_lower - mag_lower)/2 + mag_lower;
            midline1 = arcLength(midline1,'arcLen',0);                
            %%%%%%%%%%%%%%%%%        
            % get midlines
            %%%%%%%%%%%%%%%%%
            midline2 = (matched_upper - mag_upper)/2 + mag_upper;
            midline2 = arcLength(midline2,'arcLen',0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VECTOR FIELDS - tangent and normal bundle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%
            % get frames and store vector fields for midline 1
            %%%%%%%%%%%%%%%%%
            frame_bundle1 = igetFrame(midline1,5);
            %%%%%%%%%%%%%%%%%   
            % get frames and store vector fields for midline 2 and midline 1
            %%%%%%%%%%%%%%%%%
            frame_bundle2 = igetFrame(midline2,5);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OBSERVATIONS - angle and length
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%
            % extract needed
            %%%%%%%%%%%%%%%%%
            [L1 T1] = observeMidlineFibre(midline1);
            [L2 T2] = observeMidlineFibre(midline2);
            %%%%%%%%%%%%%%%%%
            % make measurements
            %%%%%%%%%%%%%%%%%
            tipAngle1 = atan2(-T1(1),-T1(2));
            tipAngle2 = atan2(-T2(1),-T2(2));    
            tipAngle = [tipAngle1 tipAngle2];
            tipAngle = mean(tipAngle);
            length = mean([L1 L2]);
            T = buildAffine(flipud(T1),fliplr(midline1(1,:)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STORE - structure
        %%%%%%%%%%%%%%%%%%
        db = cell(0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 1;
        db{IDX}.class = 'phytoACurve';
        db{IDX}.data = [fliplr(midline1) ones(size(midline1,1),1)];
        db{IDX}.cl = 'r';
        db{IDX}.icmd = 'setMidline';
        db{IDX}.name = 'midline1';
        db{IDX}.field = 'midline_1';
        db{IDX}.dims = size(midline1);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 2;
        db{IDX}.class = 'phytoACurve';
        db{IDX}.data = [fliplr(midline2) ones(size(midline2,1),1)];
        db{IDX}.cl = 'r';
        db{IDX}.icmd = 'setContour';
        db{IDX}.name = 'midline2';
        db{IDX}.field = 'midline_2';
        db{IDX}.dims = size(midline2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 7;
        db{IDX}.class = 'phytoACurve';
        db{IDX}.data = [fliplr(sec_d) ones(size(sec_d,1),1)];
        db{IDX}.cl = 'g';
        db{IDX}.icmd = [];
        db{IDX}.name = 'contour';
        db{IDX}.field = 'contour';
        db{IDX}.dims = size(contour);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 3;
        db{IDX}.class = 'phytoApoint';
        db{IDX}.data = [fliplr(tip) 1];
        db{IDX}.cl = 'g';
        db{IDX}.icmd = 'setTip';
        db{IDX}.name = 'tip_point';
        db{IDX}.field = 'tip_point';
        db{IDX}.dims = size(tip);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 4;
        db{IDX}.class = 'phytoApoint';
        db{IDX}.data = [fliplr(base) 1];
        db{IDX}.cl = 'r';
        db{IDX}.icmd = 'setBase';
        db{IDX}.name = 'base_point';
        db{IDX}.field = 'base_point';
        db{IDX}.dims = size(base);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 5;
        db{IDX}.class = 'phytoApoint';
        db{IDX}.data = [fliplr(upperBase) 1];
        db{IDX}.cl = 'r';
        db{IDX}.icmd = 'setUpperBase';
        db{IDX}.name = 'upper_base_point';
        db{IDX}.field = 'upper_base_point';
        db{IDX}.dims = size(upperBase);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 6;
        db{IDX}.class = 'phytoApoint';
        db{IDX}.data = [fliplr(lowerBase) 1];
        db{IDX}.cl = 'r';
        db{IDX}.icmd = 'setLowerBase';
        db{IDX}.name = 'lower_base_point';
        db{IDX}.field = 'lower_base_point';
        db{IDX}.dims = size(lowerBase);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IDX = 8;
        db{IDX}.class = 'phytoAaffine';
        db{IDX}.data = T;
        db{IDX}.cl = 'm';
        db{IDX}.icmd = [];
        db{IDX}.name = 'tip_reference_frame';
        db{IDX}.field = 'tip_reference_frame';
        db{IDX}.dims = size(T);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        %%% report
        if report.verbose;fprintf(['\nDone in: ' num2str(etime(clock,tm)) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    catch ME
        stop = 1;
    end
end




%{

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % STORE - return in vector bundle format - direct product
        %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        out = cell(0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store midline 1
        out{end+1} = mat2pct(midline1,'phytoCurve');
        out{end}.setHname('midline 1');
        % store midline 2
        out{end+1} = mat2pct(midline2,'phytoCurve');
        out{end}.setHname('midline 2');
        % store contour
        out{end+1}= mat2pct(sec_d,'phytoCurve');
        out{end}.setHname('contour');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store upper base
        out{end+1} = mat2pct(base,'phytoPoint');
        out{end}.setHname('upper_base');
        % store tip
        out{end+1} = mat2pct(tip,'phytoPoint');
        out{end}.setHname('tip');
        % store lower base
        out{end+1} = mat2pct(lowerBase,'phytoPoint');
        out{end}.setHname('lower_base');
        % store base point
        out{end+1} = mat2pct(base,'phytoPoint');
        out{end}.setHname('base');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store tip frame
        out{end+1} = mat2pct(T,'affine2D');
        out{end}.setHname('tip frame');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store angle measurement
        out{end+1} = mat2pct(tipAngle,'observation');
        out{end}.setHname('tip angle');
        % store angle measurement
        out{end+1} = mat2pct(length,'observation');
        out{end}.setHname('midline length');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}