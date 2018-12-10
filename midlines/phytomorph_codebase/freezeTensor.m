function [fT] = freezeTensor(T,isSet)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fT = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % place tensor in cell
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~iscell(T)
        T = {T};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % measure grade
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grade = numel(T);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store grade
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 1 | isSet == false
        isSet = false;
        setSize = 1;
    else
        for g = 1:grade
            tmpSZ = size(T{g});
            grade_check(g) = tmpSZ(end);
        end
        if ~all(grade_check == grade_check(1))
            return 
        end
        setSize = grade_check(1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store grade
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fT = repmat(grade,[1 setSize]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each rank
    for g = 1:grade
        % rank of tensor
        rnk(g) = ndims(T{g});
        % handle last rank is set index
        if isSet;rnk(g) = rnk(g) - 1;end
        % dim of tensor
        dm{g} = size(T{g});
        if isSet;dm{g}(end) = [];end
        % tmp size
        tmpSZ = size(T{g});
        newSZ = prod(tmpSZ);
        % if set then reshape different
        if isSet;newSZ = [newSZ/tmpSZ(end) tmpSZ(end)];end
        newSZ = [newSZ 1];
        % call reshape
        tmpV = reshape(T{g},newSZ);
        % make header
        header = repmat([rnk(g);dm{g}(:)],[1 size(tmpV,2)]);
        % linearize
        fT = [fT;header;tmpV];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end




%{
T = {rand(4,5),rand(3,4,5,6)};
fT = freezeTensor(T);
%}













