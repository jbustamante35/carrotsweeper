function [] = getRepresentation_ver0(f_fileList)
    for e = 1:numel(f_fileList)
        [index{e}] = selectIndex_ver0(f_fileList{e},FI,T);
        t_fileList = featureFunctionBank.generateFeatureFiles_forFile(fileList{e},oPort,keys(1:2));
        [index{e}] = intersectIndex_ver0(index{e},t_fileList);
        [data] = loadFeatures(t_fileList,index{e});
    end
end