function [k] = makeKey(k)
    warning off
    k = lower(k);
    k = strrep(k, '.','');
    k = strrep(k, '-','');
    k = strrep(k, ' ','');
    k = strrep(k, '[','');
    k = strrep(k, ']','');
    k = strrep(k, '_','');
    k = strrep(k, ')','');
    k = strrep(k, '(','');
    if ~isempty(k)
        if strcmp(k(end),'a')
            k = k(1:end-1);
        end
        if strcmp(k(end),'c')
            k = k(1:end-1);
        end
        if strcmp(k(end),'d')
            k(end) = 'b';
        end
    end

end