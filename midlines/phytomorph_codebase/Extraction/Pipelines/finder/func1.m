function [data] = func1(image,para)
    data = [];
    for k = 1:size(image,3)
        data = cat(3,data,im2colF(image(:,:,k),para.patchSize(1:2),[1 1]));
    end
    data = sum(sum(para.func(single(data)),2),1);
end