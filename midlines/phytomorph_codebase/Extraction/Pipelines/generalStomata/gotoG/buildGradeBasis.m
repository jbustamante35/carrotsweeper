function [basisU,basisE] = buildGradeBasis(gradeFileList,toLoad)

   
    [data] = myRemoteLoader(gradeFileList{1},'T');
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % compute to load and preallocate
    %%%%%%%%%%%%%%%%%%%%%%%%
    toLoad = min(toLoad,numel(gradeFileList));
    basisStack = zeros(size(data.T,1),toLoad*size(data.T,2));
       
        
    %%%%%%%%%%%%%%%%%%%%%%%%
    % init pointers
    skip = size(data.T,2);
    str = 1;
    stp = skip;
    %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start loading data for basis buildout.\n']);
    for e = 1:toLoad
        fprintf(['start loading:' num2str(e) 'data for basis buildout.\n']);
        data = myRemoteLoader(gradeFileList{e},'T');
        basisStack(:,str:stp) =  data.T;
        str = stp + 1;
        stp = str + skip - 1;
        fprintf(['end loading:' num2str(e) 'data for basis buildout.\n']);
    end
    fprintf(['end loading data for basis buildout.\n']);

    [basisStack] = transForm_extractToFeatureTensor(basisStack);
    
    for g = 1:3
        F = squeeze(basisStack(:,:,g,:));
        szF = size(F);
        F = reshape(F,[prod(szF(1:2)) szF(3)]);
        F = permute(F,[2 1]);
        [basisU{g},basisE{g}] = PCA_FIT_FULLws(F,13);
    end
    
    
end




