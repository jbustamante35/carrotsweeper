function [I] = subBackGrd(I,para)
    closeVALUE = para{1}.value;    
    resizeVALUE = para{2}.value;
    % if need to resize
    if resizeVALUE ~= 1
        sz = size(I);
        Iorg = I;
        I = imresize(I,resizeVALUE);
    end
    %%% obtain the background
    BK = getBackground(I,closeVALUE);
    % if need to resize
    if resizeVALUE ~= 1
        BK = imresize(BK,sz);
        I = Iorg;
    end
    %%% subtract the background
    I = I - BK;
end