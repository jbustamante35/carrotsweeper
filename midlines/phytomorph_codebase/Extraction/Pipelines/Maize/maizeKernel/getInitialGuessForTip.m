function [tipPoint dB] = getInitialGuessForTip(dB)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                getInitialGuessForTip.m  (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                cwtK_closed_imfilter.m, circshift.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                dB:      The information is needed. 
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        % smoothing parameters
        para{1} = [10:100];
        parfor e = 1:numel(dB)
            % get the curvature values for gaussian smoothed contours
            out = cwtK_closed_imfilter(dB{e},para);
            % get the size of the kurvature matrix
            szK = size(out.K);
            % stack for closed contour
            K = [out.K' out.K' out.K'];
            % smooth the curvature
            K = imfilter(K,fspecial('disk',21),'circular');
            % get the watersheed
            L = watershed(K);
            % find the unique watersheds
            UQ = unique(L);
            % init average depth of watershed
            uK = [];
            % for each watershed - get the average kurvature depth
            for u = 1:numel(UQ)
                % find the watersheed
                tmpW = L==UQ(u);
                % pad with zeros on top and bottom for removing watersheds            
                tmpW = [zeros(1,size(tmpW,2));tmpW;zeros(1,size(tmpW,2))];
                % clear watersheds attached to left and right
                tmpW = imclearborder(tmpW);
                % remove padding
                tmpW = tmpW(2:end-1,:);
                % find the water shed locations
                fidx = find(tmpW);
                % assign the average depth
                if ~isempty(fidx)
                    uK(u) = mean(K(fidx));
                else
                    uK(u) = inf;
                end
            end
            % find the watershed with the greatest average depth
            [~,midx] = min(uK);
            % subtract one  because watersheds start with 0
            midx = midx - 1;
            % find the deepest average watersheddB = dB(fidx);
            fidx = find(L==midx);
            % find the deepest point on the watershed
            [~,sidx] = min(K(fidx));
            % get the location
            midx = fidx(sidx);
            % convert to sub index
            [r,c] = ind2sub(size(K),midx);
            % find the location on the contour
            IDX = [1:szK(1) 1:szK(1) 1:szK(1)];
            sidx = IDX(c);            
            % get the tip point
            tipPoint(e,:) = dB{e}(sidx,:);
            % circ shift the contour for tip forward
            dB2{e} = circshift(dB{e},[-sidx 0]);
            % report done
            fprintf(['done with tip finding ' num2str(e) ':' num2str(numel(dB)) '\n'])
        end    
        dB = dB2;
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:getInitialGuessForTip.m******\n']);
    end
end