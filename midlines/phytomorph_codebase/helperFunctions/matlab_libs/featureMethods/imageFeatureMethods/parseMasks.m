function [BM] = parseMasks(B) 
    % here, anything that is connected to the border will be a seperate
    % object
    BM = [];
    % for each of the binary objects
    for im = 1:size(B,3)
        %%%%%%%%
        % move tensor into temp
        temp = B(:,:,im);
        %%%%%%%%
        % rotate and set to zero
        for ro = 1:4
            temp(:,1) = 0;
            temp = imrotate(temp,90);
        end
        %%%%%%%%
        % copy second column into first
        temp(:,1) = temp(:,2);
        %%%%%%%%
        % clear
        ctemp = imclearborder(temp);
        temp = xor(temp,ctemp);
        %%%%%%%%
        % make set of objects
        CC = bwconncomp(temp);        
        for ob = 1:CC.NumObjects
            Z = zeros([size(B,1) size(B,2)]);
            Z(CC.PixelIdxList{ob}) = 1;
            BM = cat(3,BM,Z);
        end
    end
end