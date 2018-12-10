function [outPort] = configure_outPort(outPort)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate unique key
    outPort.nDepth = length(outPort.mainBase);
    outPort.uniqueKey = genUQkey(char(stack.get(0).getFullFileName()),outPort.specialCharater,outPort.nDepth);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set the index for the matlab dat
    data.toStore.IDX = 1;
end