function [C] = applyRepresentation_ver0(fileList,representationKey,oPort,R)
    for e = 1:numel(fileList)
        % construct the feature file name for the feature which needs a
        % representation built
        t_fileList = featureFunctionBank.generateFeatureFiles_forFile(fileList{e},oPort,{representationKey});
        % load feature data
        fprintf(['Starting load of features \n']);tic
        [data{e} sz] = loadFeatures(t_fileList);
        fprintf(['Ending load of features:' num2str(toc) '\n']);
    end
    data = cell2mat(data);
    data = mbindVec(data);
    C = PCA_REPROJ_T(data,R.E,R.U);
    for e = 1:100
        K = reshape(C(e:(e+2),:)',[sz{1}-2*16 3]);
        imshow(K,[]);
        drawnow
    end
end