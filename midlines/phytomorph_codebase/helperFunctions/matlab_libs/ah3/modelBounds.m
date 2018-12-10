function [b] = modelBounds(data)
    uMIN = min(data,[],1);
    uMAX = max(data,[],1);            
    SIGMA = uMAX - uMIN;
    sMIN = zeros(size(SIGMA));
    %sMIN = SIGMA/10;
    sMAX = SIGMA;
    b = [[uMIN;uMAX],[sMIN;sMAX]];
end