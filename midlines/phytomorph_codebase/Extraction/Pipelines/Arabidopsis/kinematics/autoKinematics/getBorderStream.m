function [stream] = getBorderStream(I,SZ)
    stream = [];
    for r = 1:4
        stream = [stream;I(:,1:SZ)];
        %I(:,1:SZ) = [];
        I = imrotate(I,-90);
    end
    stream = im2col(stream,[SZ SZ],'sliding');
end