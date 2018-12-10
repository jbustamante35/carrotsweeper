function [y] = myOp1(m,v)
    %y = mtimesx(v',permute(m.[2 1 3]),'speed');
    for k = 1:size(m,3)
        y(:,k) = mtimesx(single(v(k,:)),single(m(:,:,k)),'speed');
    end
end