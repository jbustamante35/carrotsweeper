function [out] = cwtK_closed_peaks(J,para,dz,thresh)
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
        pK = imdilate(abs(K),ones(1,dz))==abs(K);
        tK = abs(K) > thresh;
        pIDX = find(tK.*pK);
        K = K(:,(sz(1)+1):(sz(1)+sz(1)));
        k = (pIDX > sz(1)) & (pIDX < 2*sz(1));
        pIDX = pIDX(k) - sz(1) + 1;
        % outPort
        out.K = K;
        out.pIDX = pIDX;
    catch ME
        ME.message
    end
end