function [outPort] = generate_outPort(outPort,para)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % setup outPort
    % generate constants
    outPort.basePath = para.basePath;
    outPort.nDepth = para.nDepth;
    outPort.specialCharater = '--';    
    if para.request_uniqueKey
        outPort.uniqueKey = genUQkey(para.localFileBase,outPort.specialCharater,outPort.nDepth);
    end
    
    % disk I/O
    outPort.disk = 1;                                   % true for writing to disk
    outPort.toDisk = [];                                % data to write to disk
    outPort.nDepth = 1;                                 % depth up the file tree for generating unique name
    
    % store I/O
    outPort.oStore = 0;                                 % true for writing to object store        
    outPort.toStore = @(x1,x2)objectSpool_mm(x1,x2);    % data to write to disk
    
    % zip I/O
    outPort.zip = 0;                                    % true for zipping data for pull back
end