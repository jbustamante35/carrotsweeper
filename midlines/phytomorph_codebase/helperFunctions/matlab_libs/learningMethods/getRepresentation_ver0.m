function [R] = getRepresentation_ver0(fileList,FI,T,selectKey,representationKey,vecCompSelect,oPort)
    [f_fileList] = featureFunctionBank.generateFeatureFiles_forKey(fileList,oPort,selectKey);
    parfor e = 1:numel(f_fileList)
        % construct the feature file name for the feature which needs a
        % representation built
        t_fileList = featureFunctionBank.generateFeatureFiles_forFile(fileList{e},oPort,{representationKey});
        % find the index from the select feature function
        [index] = selectIndex_ver0(f_fileList{e},FI,T);
        % find the locations which can be selected
        [index] = intersectIndex_ver0(index,t_fileList);
        fprintf(['Starting load of features \n']);tic
        % load feature data
        [data{e}] = loadFeatures(t_fileList,index);
        data{e} = data{e}(:,:,vecCompSelect);
        fprintf(['Ending load of features:' num2str(toc) '\n']);
    end
    data = cell2mat(data);
    ndata = mbindVec(data);
    [S C U E L ERR LAM] = PCA_FIT_FULL_T(ndata,size(data,1));
    R.E = E;
    R.U = U;
end