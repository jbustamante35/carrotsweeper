function [foo] = extractFunction_0(varargin)
    fprintf(['******************************************** \n']);
    fprintf(['Start function: import \n']);
    fprintf(['******************************************** \n']);
    % goal: transfer image type into featureMap
    % import data from image into 1x1 pixel paches
    % data is stored as row vector in raster format
    %[pth,nm,ext] = fileparts(varargin{1}.fileName);
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
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% preform the operation(s)
    patchSize = 0;
    data = myReader(varargin{1});
    data = data(:,:,1);
    
    SZ = size(data);
    index = 1:(SZ(1)*SZ(2));
    index = reshape(index,SZ(1:2));
    %% return feature map object
    foo = fmo();
    %foo.setHeadNodeName(nm);
    foo.setData(data(:)');
    foo.setIndexMap(index(:)');
    foo.setPatchSize(patchSize);
    foo.setRawSize(SZ);
    foo = {foo};
    fprintf(['******************************************** \n']);
    fprintf(['End function: import \n']);
    fprintf(['******************************************** \n']);
end