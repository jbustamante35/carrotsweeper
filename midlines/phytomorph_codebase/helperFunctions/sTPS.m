function [] = sTPS(D,n)

    IDX = zeros(numel(n),2);
    str = 1;
    for e = 1:numel(n)
        b = nchoosek(size(D,2),n(e));
        stp = str + b - 1;
        IDX(e,:) = [str stp];
        str = stp + 1;
    end
   
    R = zeros(size(D,1),IDX(end,2));
    for e = 1:numel(n)
        R(:,IDX(e,1):IDX(e,2)) = sTP(D,n(e));
    end

end

%{

    sTPS(rand(13,5),[1 2])

%}