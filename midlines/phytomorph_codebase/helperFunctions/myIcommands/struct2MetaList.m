function [nm] = struct2MetaList(S)
    nm = [];
    for i = 1:numel(S)
        nm = [nm '{[' S(i).key '].[' S(i).value ']}'];
    end
end