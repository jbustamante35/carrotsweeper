function [foo] = extractFunction_2(varargin)
    fprintf(['******************************************** \n']);
    fprintf(['Start function: kurvature \n']);
    fprintf(['******************************************** \n']);
    %% assign the varargin
    % 1: the data object to compute on
    % 2: the scale of kurvature operation
    % 3: the patch size post kurvature compute
    %% parpare the varargin{1};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare the varargin{1}
    if isa(varargin{1},'org.openrdf.model.impl.URIImpl')
        varargin{1} = char(varargin{1}.toString());
    end
    if ischar(varargin{1})
        fidx = strfind(varargin{1},':');
        if ~isempty(fidx)
            varargin{1} = varargin{1}((fidx(1)+2):end);
        end
        tmp = load(varargin{1},'obj');
        varargin{1} = tmp.obj;
        clear tmp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = varargin{1}.getData();
    kScale = varargin{2};
    patchSize = varargin{3};
    % get the raw size of the feature
    SZ = varargin{1}.getRawSize();
    % reshape for computation
    data = reshape(data,SZ);
    % create new index map
    indexMapFromInput = varargin{1}.getIndexMap();
    indexMapFromInput = reshape(indexMapFromInput,SZ);
    index = trimGray(indexMapFromInput,patchSize);
    newSZ = size(index);
    %% perform kurvature calculation
    para.scales.value = kScale;
    para.resize.value = 1;
    lK = surKur(data,para);
    SZ = size(lK);
    %% assign the outputs
    for e = 1:SZ(end)
        foo{e} = fmo();
        %foo{e}.setHeadNodeName(varargin{1}.getHeadNodeName());
        foo{e}.setData(im2colF(squeeze(lK(:,:,e)),2*patchSize+1,[1 1]));
        foo{e}.setIndexMap(index(:)');
        foo{e}.setPatchSize(patchSize);
        foo{e}.setRawSize(newSZ);
    end
    fprintf(['******************************************** \n']);
    fprintf(['End function: kurvature \n']);
    fprintf(['******************************************** \n']);
end