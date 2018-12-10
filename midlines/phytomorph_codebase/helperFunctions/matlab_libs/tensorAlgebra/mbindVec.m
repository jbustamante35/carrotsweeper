function [vec] = mbindVec(vec)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clip-normalization of the vector        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           vec     := vector to operate on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           vec     := resulting vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    vec = bsxfun(@minus,vec,min(vec,[],1));
    vec = bsxfun(@times,vec,max(vec,[],1).^-1);
end
