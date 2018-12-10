function [out] = cwtK_closed_imfilter(J,para)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                cwtK_closed_imfilter.m  (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                fspecial.m, gradient.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                J:      The information is needed. 
                para:      The information is needed.
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    % MUST REMOVE THE BASE
    % FIND AND SUPPRES    
    try
        for e = 1:numel(para{1})
            h0 = fspecial('gaussian',[1 5*para{1}(e)],para{1}(e));
            h0 = h0 / sum(h0);
            % calculate curvature            
            tmp = imfilter(J,h0','circular');
            d1X = gradient(tmp')';
            d2X = gradient(gradient(tmp'))';            
            K(:,e) = (d1X(:,1).*d2X(:,2) - d1X(:,2).*d2X(:,1)).*(d1X(:,1).^2 + d1X(:,2).^2).^-3/2;    
        end
        
        % outPort
        out.K = K;

    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:cwtK_closed_imfilter.m******\n']);
    end

end