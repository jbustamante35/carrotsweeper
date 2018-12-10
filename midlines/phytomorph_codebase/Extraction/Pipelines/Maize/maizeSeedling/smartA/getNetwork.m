function [funcObject] = getNetwork(W,numC,netS,Y)
    % W: size of [x y 3 trials
    % get size
    %sz = size(W);
    if ndims(W) ~= 2
        % permute to make choice about columns
        W = permute(W,[1 3 2 4]);
        sz = size(W);
        % stack color panels 
        W = reshape(W,[prod(sz(1:2)) prod(sz(3:4))]);
    end
    % pca fit
    [wU wE] = PCA_FIT_FULL_Tws(W,numC);
    wC = PCA_REPROJ_T(W,wE,wU);
    % clear W
    clear W
    %% make pattern net
    net = patternnet(netS);
    net = train(net,wC,logical(Y),'useParallel','yes');
    funcObject.net = net;
    funcObject.U = wU;
    funcObject.E = wE;
end