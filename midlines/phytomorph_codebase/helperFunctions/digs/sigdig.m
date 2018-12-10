function [dirList] = sigdig(P,type)
    [dirList] = mydir(P);
    dirList =  struct2list(dirList,'name');
    idx = isType(dirList,type);
    dirList = dirList(idx);
end