function [out] = cwtK_closed(J,para)
    % MUST REMOVE THE BASE
    % FIND AND SUPPRES    
    try
        sz = size(J);
        J = [J;J;J];
        % calculate curvature
        d1X1 = cwt(J(:,1),para{1},'gaus1');
        d1X2 = cwt(J(:,2),para{1},'gaus1');
        d2X1 = cwt(J(:,1),para{1},'gaus2');
        d2X2 = cwt(J(:,2),para{1},'gaus2');
        K = (d1X1.*d2X2 - d1X2.*d2X1).*(d1X1.^2 + d1X2.^2).^-3/2;
        K = K(:,sz(1)+1:sz(1)+sz(1));
        
        % outPort
        out.K = K;

    catch ME
        ME.message
    end
end