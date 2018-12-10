function [N] = imgNormal(I,para)
    gradPara = para{4};                            % = gradient var
    D = diffmethod(I,gradPara);                     % get the information for the structure tensor
    N = cat(3,-D(:,:,1),-D(:,:,2),ones(size(I)));
end