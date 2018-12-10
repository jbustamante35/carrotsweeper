function [ret] = tensorEdge(T,para)
    %%%%%%%%%%%%%%%%%%
    % min value
    mH = min(T(:,para{1}),[],1);
    % select those which match this value
    fidx = T(:,para{1}) == mH;
    % return the points
    ret.T = T;
    ret.fidx = fidx;
    %%%%%%%%%%%%%%%%%%
end