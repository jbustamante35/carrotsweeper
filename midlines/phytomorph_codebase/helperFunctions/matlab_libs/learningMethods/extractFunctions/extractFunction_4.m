function [foo] = extractFunction_4(varargin)
    %% assign the varargin
    % 1: the data object to compute on
    % 2: the basis to compute the error with respect to
    % 3: the patch size of the error
    data = varargin{1}.getData();
    basisInformation = varargin{2};
    patchSize = varargin{3};
    ncomps = 20;
    % the size of the old feature map
    SZ = varargin{1}.getRawSize();
    % create new index map
    indexMapFromInput = varargin{1}.getIndexMap();
    indexMapFromInput = reshape(indexMapFromInput,SZ);
    index = trimGray(indexMapFromInput,patchSize);
    newSZ = size(index);
    % calculate feature map    
    [C ERR] = PCA_REPROJ_T(data,basisInformation.E(:,ncomps),basisInformation.U);
    ERR = reshape(ERR,SZ);
    %% assign the outputs
    for e = 1:1
        foo{e} = fmo();
        foo{e}.setHeadNodeName(varargin{1}.getHeadNodeName());
        foo{e}.setData(im2colF(ERR,2*patchSize+1,[1 1]));
        foo{e}.setIndexMap(index(:)');
        foo{e}.setPatchSize(patchSize);
        foo{e}.setRawSize(newSZ);
    end
end