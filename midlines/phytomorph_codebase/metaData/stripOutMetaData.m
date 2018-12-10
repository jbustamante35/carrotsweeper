function [md] = stripOutMetaData(string)
    fidx_str = strfind(string,'{');
    fidx_end = strfind(string,'}');
    for e = 1:numel(fidx_str)
        sub = string((fidx_str(e)+1):(fidx_end(e)-1));
        strD = strfind(sub,'_');
        key = sub(1:(strD(1)-1));
        value = sub((strD(1)+1):end);
        md.(key)=  value;
    end
end