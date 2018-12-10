function [] = setMetaGfromStruct(L,S)
    setMetafromStruct(L{1},S);
    for i = 2:size(L,1)
        cmd = ['imeta cp -d -d ' L{1} ' ' L{i}];
        [status query] = system(cmd);
    end
end