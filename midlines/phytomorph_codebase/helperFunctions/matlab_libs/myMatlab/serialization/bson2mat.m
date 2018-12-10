function [out] = bson2mat(in)
    out = masterSwitch(in);
end
%%%%%%
% recursive call to this
function [out] = masterSwitch(in)
    switch in.get('class')
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
    for e = 1:in.get('data').size()
        out{e} = masterSwitch(in.get('data').get(num2str(e)));        
    end
end        
%%%%%%
% convert struct
function [out] = convertStruct(in)    
    keySet = in.get('data').keySet();
    itr = keySet.iterator();
    while itr.hasNext()        
        key = itr.next();
        out.(char(key)) = masterSwitch(in.get('data').get(key));
    end
end
%%%%%%
% convert a matlab charater to String
function [out] = convertChar(in)
    out = char(in.get('data').get('value'));
end
%%%%%%
% convert tensor to
function [out] = convertTensor(in)    
    sz = in.get('size');
    out = zeros(prod(sz));
    data = in.get('data'); 
    out = cell2mat(data.toArray.cell);
    out = reshape(out,sz');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

    T{1}.what = 'hello';
    T{2} = 4;
    bT = mat2bson(T);
    mT = bson2mat(bT);
%}