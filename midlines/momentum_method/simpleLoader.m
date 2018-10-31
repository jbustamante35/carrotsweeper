function [DATA] = simpleLoader(searchStr,FileList,FR)
    ret = find(strcontains(FileList,searchStr));
    DATA = [];
    for e = 1:numel(ret)
        fn = FileList{ret(e)};
        D = csvread(fn);
        if size(D,1) >= FR
            DATA = [DATA D(1:FR,:)];
        end
    end
end