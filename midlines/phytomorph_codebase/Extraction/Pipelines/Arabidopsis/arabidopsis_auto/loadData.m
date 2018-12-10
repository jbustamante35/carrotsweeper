function [out] = loadData(file,NP,SNIP,NPK)
    % load the mat file    
    out = load(file,'out');
    % assign to the proper variable
    out.frameData = out.out;
    out = rmfield(out,'out');
    % measure angle,length,kurvature,growthrate
    out = measureStruct(out,NP,NPK,SNIP,[301 1]);
    out.growthRate = gradient(out.length')';
    [pth nm ext] = fileparts(file);
    % sort fields
    out = orderfields(out);
    % attach data ID    
    out.dataID = nm;
end