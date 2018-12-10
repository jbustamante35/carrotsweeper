function [inStruct] = myQRparse(stringIN,inStruct)   
    stringIN = stringIN{1};
    fidx1 = strfind(stringIN,'{');
    fidx2 = strfind(stringIN,'}');
    if size(inStruct,1) == 0
        init = true;
    else
        init = false;
    end
    for e = 1:numel(fidx1)
        datum = stringIN((fidx1(e)+1):(fidx2(e)-1));
        mysep = strfind(datum,'_');
        key = datum(1:(mysep(1)-1));
        value = datum(mysep(1)+1:end);
        
        if init
            try
                datenum(value);
            catch
                if ~isempty(str2num(value))
                    value = str2num(value);
                end
            end
        else
            if ~iscell(inStruct.(key))
                value = str2num(value);
            end
        end
        
       
        
        if init
            if ischar(value)
                inStruct.(key){1,1} = value;
            else
                inStruct.(key)(1,1) = value;
            end
        else
            if ischar(value)
                 inStruct.(key){end+1,1} = value;
            else
                inStruct.(key)(end+1,1) = value;
            end
        end
    end
end