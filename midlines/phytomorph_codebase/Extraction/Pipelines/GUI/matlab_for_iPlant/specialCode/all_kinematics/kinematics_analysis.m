function [profile] = kinematics_analysis(rts,para)
    try
        switch para.kinematicsType
            case 'steadyState'
                [profile] = rts.steadyStateFlowAnalysis();                
            case 'nonSteadyState'
                [profile] = rts.steadyStateFlowAnalysis();
        end
    catch ME
        ME;
    end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% strain analysis : sliding window and expectation values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uY uX] = strainAnalysis_type0(X,Y,Xi,window,disp)
    % obtain the unqiue domain    
    [X,sidx] = unique(X);
    % obtain the corresponding coDomain
    Y = Y(sidx);
    % sort domain and doDomain
    [X sidx] = sort(X);
    Y = Y(sidx);
    % for each 
    for e = 1:numel(Xi)
        % find the points in the domain
        fidx = find((X > (Xi(e) - window)) &  (X < (Xi(e) + window)));    
        EM(e) = isempty(fidx);
        % obtain the corresponding coDomain points
        iY = Y(fidx);
        % obtain the expected value (funny joke, mean or average)
        uY(e) = mean(iY);
        uX(e) = mean(X(fidx));
        % disp
        if disp
            plot(X,Y,'b.');
            hold on;
            plot(X(fidx),Y(fidx),'r.');
            plot(Xi(1:e),iY(1:e),'g');
            hold off
            drawnow
        end
    end
end
