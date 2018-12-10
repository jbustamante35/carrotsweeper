function [foo] = extractFunction_1(varargin)
    fprintf(['******************************************** \n']);
    fprintf(['Start function: sliding window \n']);
    fprintf(['******************************************** \n']);
    % sliding window operation from raw image data
    % make sliding window operation as a feature
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
    %% run the operation
    % set values from varargin
    data = varargin{1}.getData();
    patchSize = varargin{2};
    SZ = varargin{1}.getRawSize();
    data = reshape(data,SZ);
    
    % create new index map
    indexMapFromInput = varargin{1}.getIndexMap();
    indexMapFromInput = reshape(indexMapFromInput,SZ);
    index = trimGray(indexMapFromInput,patchSize);
    newSZ = size(index);
    %% return feature map object
    foo = fmo();
    %foo.setHeadNodeName(varargin{1}.getHeadNodeName());
    foo.setData(im2colF(data,2*patchSize+1,[1 1]));
    foo.setIndexMap(index(:)');
    foo.setPatchSize(patchSize);
    foo.setRawSize(newSZ);
    foo = {foo};
    fprintf(['******************************************** \n']);
    fprintf(['End function: sliding window \n']);
    fprintf(['******************************************** \n']);
end