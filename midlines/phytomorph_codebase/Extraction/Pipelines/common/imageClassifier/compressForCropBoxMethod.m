function [toCompress] = compressForCropBoxMethod(toCompress,U,E)
    for e = 1:numel(toCompress)
        toCompress{e} = PCA_REPROJ_T(toCompress{e},E,U);
    end
end