function [ukey] = generateUniqueKey()
    ukey = strrep([num2str(now) num2str(rand(1,1))],'.','');
end