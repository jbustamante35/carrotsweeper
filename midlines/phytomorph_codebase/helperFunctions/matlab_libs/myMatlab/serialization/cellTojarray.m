function [out] = cellTojarray(varargin)
    % import libs
    import phytoG.locked.BdataStructures.implementations.directMongo.nan.*;
    %%%%%%%%%%%%%%%%%%%%%%
    % array type
    if nargin == 1
        type = 'java.lang.Object';
    elseif nargin == 2
        type = varargin{2};
    end
    %%%%%%%%%%%%%%%%%%%%%%
    % return object
    out = javaArray(type,numel(varargin{1}));
    %%%%%%%%%%%%%%%%%%%%%%
    % put into return object
    for e = 1:numel(varargin{1})
        cmd = ['out(e) = ' type '(varargin{1}{e});'];
        eval(cmd);
    end
    %%%%%%%%%%%%%%%%%%%%%%
    % convert to B2DBL
    out = B2DBL(out);
end