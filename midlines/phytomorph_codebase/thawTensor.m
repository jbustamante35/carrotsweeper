function [T] = thawTensor(fT,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get pointer list
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ptrList,dm] = calculatePTRs(fT(:,1));
    % if nargin is one, then thaw all
    if nargin == 1
        n = 1:size(ptrList,1);
    end
    % for each grade
    for g = 1:numel(n)
        % get the first tensor for preallocation
        V = fT(ptrList(1,n(g)):ptrList(2,n(g)),1);
        % for each tensor
        V = fT(ptrList(1,n(g)):ptrList(2,n(g)),:);
        % restore size
        V = reshape(V,[dm{n(g)}' size(V,2)]);
        T{g} = V;
    end
    if numel(n) == 1
        T = T{1};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

function [ptrList,dm] = calculatePTRs(fT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % measure grade
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ptrList = [];
    grade = fT(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each rank
    ptr = 2;
    for g = 1:grade
        % get rank
        rnk = fT(ptr);
        % get dim of tensor
        dm{g} = fT(ptr+1:ptr+rnk);
        % set str
        str = ptr+rnk+1;
        % set stp
        stp = ptr+rnk+prod(dm{g});
        % stack pointer list
        ptrList = [ptrList [str;stp]];
        % increment
        ptr = stp + 1;
    end
end
%{
T = {rand(4,5),rand(3,4,5,6)};
fT = freezeTensor(T);
rT = thawTensor(fT);
rT{1} == T{1}
rT{2} == T{2}


T1 = {rand(4,5),rand(3,4,5,6)};
T2 = {rand(4,5),rand(3,4,5,6)};
fT1 = freezeTensor(T1);
fT2 = freezeTensor(T2);
rTW = thawTensor([fT1;fT2]);
rTW{1} == T1{1}
rTW{2} == T1{2}





%}













