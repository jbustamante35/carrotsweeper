function [tdB] = findTipPoints_fast(dB,B,I)
    % all warnings are not to be shown.
    warning off
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                findTipPoints_fast.m  (Inputs are relative to 1200dpi)
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
        tdB = dB;
        SEGSIZE = 601;
        [E,U] = generateCurveSegments(dB,SEGSIZE,10);

        disp = 0;
        %%%%%%%%%%%%%%%%%%%%%%%
        % compile and come back to original path
        % this line added to compile it properly.-Lee
        %%%%%%%%%%%%%%%%%%%%%%%
        parfor k = 1:numel(dB)
            fprintf(['starting creating measure tensor for:' num2str(k) ':' num2str(numel(dB)) '\n']);
            tm = clock;
            M{k} = measureSingleContourMetrics_fast(dB{k},B,E,U,SEGSIZE,10);
            fprintf(['ending creating measure tensor for:' num2str(k) ':' num2str(numel(dB)) ' in ' num2str(etime(clock,tm)) '\n']);
        end
        shiftAmount = zeros(1,numel(tdB));
        rep = 1;
        while ~(rep ~= 1 && all(shiftAmount(end,:) == 0))
            tM = [];

            for k = 1:numel(dB)
                tM(k,:) = M{k}(1,:);
            end


            H = {};
            for e = 1:size(tM,2)
                H{e} = getColumnDistribution(tM(:,e));
            end



            %imshow(I,[]);
            %hold on
            if disp
                imshow(I,[]);
                hold on;
            end
            for e = 1:numel(dB)
                MLE = lookUpLOG(H,M{e});
                MLE = sum(MLE,2);
                [J,nidx] = max(MLE);
                toShift = -(nidx-1);
                shiftAmount(rep,e) = toShift;
                if disp
                    plot(tdB{e}(:,2),tdB{e}(:,1),'r');
                    hold on
                    plot(tdB{e}(1,2),tdB{e}(1,1),'c*');
                    plot(tdB{e}(nidx,2),tdB{e}(nidx,1),'b*');
                    drawnow
                end
                tdB{e} = circshift(tdB{e},[toShift 0]);
                M{e} = circshift(M{e},[toShift 0]);
            end
            rep = rep + 1;
        end
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:findTipPoints_fast.m******\n']);
    end
end
