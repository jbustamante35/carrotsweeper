function [C] = diskCOV(fileList,U,toVecFunc,selVecFunc)
    z = 0;
    % for vector source
    for e = 1:numel(fileList)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load draw(s) from source
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ischar(fileList{e})
            I = imread(fileList{e});
        else
            I = fileList{e};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % apply func to map draw to vector(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % map to vectors
        I = toVecFunc(I,e);
        % select vectors
        I = selVecFunc(I,e);
        % init if first one
        if e == 1
            C = zeros(size(I,1),size(I,1));
        end
        % subtract mean
        I = bsxfun(@minus,I,U);
        % create COV
        C = C + I*I';
        % new image size
        z = z + size(I,2);
    end
    C = C * z.^-1;
end
