function [D] = cornerstrength(D)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grade the corner
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           D   := this is a mixture of the information in D to grade
    %                  the corner    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           D   := the mixture for corner is stacked onto D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D = cat(3,D,(D(:,:,3).*D(:,:,4) - D(:,:,5).^2).*(D(:,:,3) + D(:,:,4) + eps).^-1);
end
