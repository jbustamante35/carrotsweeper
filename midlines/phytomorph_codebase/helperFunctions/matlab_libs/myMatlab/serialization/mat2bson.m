function [out] = mat2bson(in)
    out = masterSwitch(in);
end
%%%%%%
% recursive call to this
function [out] = masterSwitch(in)
    switch class(in)
        case 'char'
            out = convertChar(in);
        case 'double'
            out = convertTensor(in);
        case 'uint8'
            out = convertTensor(in);
        case 'single'
            out = convertTensor(in);
        case 'cell'
            out = convertCell(in);
        case 'struct'
            out = convertStruct(in);
    end
end
%%%%%%
% convert cell
function [out] = convertCell(in)
    import com.mongodb.*;
    out = BasicDBObject();
    data = BasicDBObject();
    for e = 1:numel(in)
        data.append(num2str(e),masterSwitch(in{e}));
    end    
    out.append('class',class(in));
    out.append('data',data);
end        
%%%%%%
% convert struct
function [out] = convertStruct(in)
    import com.mongodb.*;
    out = BasicDBObject();
    data = BasicDBObject();
    f = fieldnames(in);
    for e = 1:numel(f)
        data.append(f{e},masterSwitch(in.(f{e})));
    end
    out.append('class',class(in));
    out.append('data',data);
end
%%%%%%
% convert a matlab charater to String
function [out] = convertChar(in)
    import com.mongodb.*;
    import java.lang.*;
    out = BasicDBObject();
    data = BasicDBObject();
    data.append('value',java.lang.String(in));    
    out.append('class',class(in));
    out.append('data',data);
end
%%%%%%
% convert tensor to 
function [out] = convertTensor(in)
    import com.mongodb.*;
    import java.lang.*;
    import java.util.*;
    % size
    SZ = size(in);
    % data
    data = Vector(prod(SZ));
    for e = 1:numel(in)
        %data(e) = java.lang.Double(in(e));
        data.add(in(e));
    end
    out = BasicDBObject();
    out.append('data',data);
    out.append('class',class(in));
    out.append('size',SZ);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

    T{1}.what = 'hello';
    T{2} = 4;
    bT = mat2bson(T);

    R = rand(10000,1);
    bR = mat2bson(R);
    mT = bson2mat(bR);
%}