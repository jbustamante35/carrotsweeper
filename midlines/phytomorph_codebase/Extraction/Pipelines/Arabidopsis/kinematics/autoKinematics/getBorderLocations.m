function [stream] = getBorderLocations(I,SZ)
    I = reshape(1:numel(I),size(I));
    stream = [];
    for r = 1:4
        stream = [stream;I(:,1:SZ)];
        %I(:,1:SZ) = [];
        I = imrotate(I,-90);
    end
    HW = (SZ-1)/2+1;
    stream = stream(HW:end-HW+1,(SZ-1)/2);
end