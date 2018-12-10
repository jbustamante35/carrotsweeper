function [E] = structExtract(S)
    F = fieldnames(S);
    for i = 1:size(F,1)
        cmd = [F{i} '= S.(F{i})'];
        eval(cmd);
    end
end
