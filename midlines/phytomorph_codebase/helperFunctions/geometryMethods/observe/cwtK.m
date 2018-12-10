function [out] = cwtK(J,para)
    % MUST REMOVE THE BASE
    % FIND AND SUPPRES    
    try
        % calculate curvature
        d1X1 = cwt(J(:,1),para{1},'gaus1');
        d1X2 = cwt(J(:,2),para{1},'gaus1');
        d2X1 = cwt(J(:,1),para{1},'gaus2');
        d2X2 = cwt(J(:,2),para{1},'gaus2');
        K = (d1X1.*d2X2 - d1X2.*d2X1).*(d1X1.^2 + d1X2.^2).^-3/2;
        K(:,1:3*max(para{1})) = 0;
        K(:,end-3*max(para{1}):end) = 0;
        
        % outPort
        out.K = K;
        out.baseSize = sum(J(:,2)==1);
        out.J = J;
    catch ME
        ME.message;
    end
end