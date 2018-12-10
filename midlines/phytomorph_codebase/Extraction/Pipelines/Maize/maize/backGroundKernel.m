function [uB] = backGroundKernel(I)
    uB = -1;
    try
        backGround = I(1:100,:);
        uB = mean(backGround(:));
    catch
        
    end
end