function [c] = tensorCombine_0(T)
    %%%%%%%%%%
    % combine
    %c = abs(T(:,:,1)-T(:,:,2));
    %c = (T(:,:,1) - T(:,:,2)).*((10^14)*T(:,:,1).*T(:,:,2)).^-1;
    c = T(:,:,2);
    c = bindVec(c);
    % normalize
    %c = bindVec(c);
    %%%%%%%%%%
    % 
end