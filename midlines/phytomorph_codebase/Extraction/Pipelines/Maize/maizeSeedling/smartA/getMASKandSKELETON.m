function [connectedMASK,MASK,SKELETON] = getMASKandSKELETON(I)
    for e = 1:numel(I)
        tic;
        fprintf(['starting:: get mask \n']);
        % get the mask
        MASK{e} = getMASK_ver0(I{e});
        % connect the mask
        connectedMASK{e} = connectPlant(MASK{e});
        % get the largest object in the mask
        connectedMASK{e} = bwlarge(connectedMASK{e});
        fprintf(['ending:: get mask::' num2str(toc) '\n']);
        tic
        fprintf(['starting:: get skeleton \n']);
        SKELETON{e} = getPlantSkelton(MASK{e});
        fprintf(['ending:: get skeleton::' num2str(toc) '\n']);
    end
end