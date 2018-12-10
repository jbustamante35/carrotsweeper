function [tmpD] = cropAndtrim(I,R,topTRIM)
    for e = 1:numel(R)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop a vertical strip for each container
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting: crop the strip\n']);
        % crop the strip
        tmpD{e} = imcrop(I,round(R(e).BoundingBox));
        % next boundingbox
        nR(e).BoundingBox = R(e).BoundingBox;
        % trim the top
        tmpD{e}(1:topTRIM:end,:,:) = [];
        fprintf(['ending: crop the strip\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop a vertical strip for each container
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end