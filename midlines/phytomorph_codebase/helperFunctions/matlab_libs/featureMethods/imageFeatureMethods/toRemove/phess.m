function [N] = phess(I,para)
    gradPara = para{4};
    D = diffmethod(I,gradPara);
    A = diffmethod(D(:,:,1),gradPara);
    B = diffmethod(D(:,:,2),gradPara);
    N = cat(3,A(:,:,1),.5*(A(:,:,2) + B(:,:,1)),.5*(A(:,:,2) + B(:,:,1)),B(:,:,2));
end
