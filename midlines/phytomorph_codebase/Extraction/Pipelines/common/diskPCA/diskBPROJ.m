function [C] = diskBPROJ(fileList,U,E,toVecFunc,selVecFunc)
    for e = 1:numel(fileList)
        if ischar(fileList{e})
            I = double(imread(fileList{e}));
        else
            I = fileList{e};
        end
        % apply to source load to get vector(s)
        I = toVecFunc(I,e);
        % select vectors
        I = selVecFunc(I,e);
        % subtract and project
        C(:,:,e) = (bsxfun(@minus,I,U)'*E)';
    end
    % squeeze
    C = squeeze(C);
end