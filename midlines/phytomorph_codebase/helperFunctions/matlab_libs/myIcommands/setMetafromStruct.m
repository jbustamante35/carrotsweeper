function [] = setMetafromStruct(L,S)
    F = fieldnames(S);
    for i = 1:size(F,1)
        cmd = ['imeta add -d ' L ' ' F{i} ' ' S.(F{i}) ' '];
        [status query] = system(cmd);
        %cmd = ['imeta ls -d ' L ];
    end
end