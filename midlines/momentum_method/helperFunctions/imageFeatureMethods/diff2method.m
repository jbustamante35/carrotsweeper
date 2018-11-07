function [Do] = diff2method(D,para)
    Do = zeros(size(D,1),size(D,2),4);    
    for d = 1:2
        st = (d-1)*2+1;
        ed = st+1;
        Do(:,:,st:ed) = diffmethod(D(:,:,d),para);
    end
end