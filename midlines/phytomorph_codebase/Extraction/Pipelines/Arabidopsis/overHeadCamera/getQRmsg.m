function [msg] = getQRmsg(CI,LI,V)
    for e = 1:numel(LI)
        label = zeros(size(LI{e}));
        for l = 1:numel(V{4})
            label = label | bwlarge(LI{e} == V{4}(l));
        end
        R = regionprops(label,'BoundingBox');
        d = imcrop(CI{e},R(1).BoundingBox);
        msg{e} = decode_qr(d);
    end
end