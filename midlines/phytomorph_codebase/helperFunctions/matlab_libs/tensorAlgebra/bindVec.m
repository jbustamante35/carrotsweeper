function [vec] = bindVec(vec)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clip-normalization of the vector        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           vec     := vector to operate on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           vec     := resulting vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    vec = vec - min(vec(:));
    vec = vec / max(vec(:));
end
