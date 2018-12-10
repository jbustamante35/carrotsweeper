function [foo] = extractFunction_3(varargin)
    fprintf(['******************************************** \n']);
    fprintf(['Start function: generate target \n']);
    fprintf(['******************************************** \n']);
    %% parpare the varargin{1};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare the varargin{1}
    toPrepare = 2;
    for e = 1:toPrepare
        if isa(varargin{e},'org.openrdf.model.impl.URIImpl')
            varargin{e} = char(varargin{e}.toString());
        end
        if ischar(varargin{e})
            fidx = strfind(varargin{e},':');
            if ~isempty(fidx)
                varargin{e} = varargin{e}((fidx(1)+2):end);
            end
            tmp = load(varargin{e},'obj');
            varargin{e} = tmp.obj;
            clear tmp;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dilateSZ = varargin{3};
    SZ = varargin{1}.getRawSize();
    %%
    Z = zeros(SZ);
    for p = 1:size(varargin{2}.pointSet,1)
       Z(varargin{2}.pointSet(p,1),varargin{2}.pointSet(p,2)) = 1;
    end
    %%
    Z = imdilate(Z,strel('disk',dilateSZ,0));
    %% index
    index = 1:(SZ(1)*SZ(2));
    index = reshape(index,SZ(1:2));
     %% return feature map object
    foo = fmo();
    %foo.setHeadNodeName(varargin{1}.getHeadNodeName());
    foo.setData(Z(:)');
    foo.setIndexMap(index(:)');
    foo.setPatchSize(0);
    foo.setRawSize(SZ);
    foo = {foo};
    fprintf(['******************************************** \n']);
    fprintf(['End function: generate target \n']);
    fprintf(['******************************************** \n']);
end