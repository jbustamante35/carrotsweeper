function [U,C] = diskMEANandCOV(fileList,toVecFunc,selVecFunc)
    % 
    U = diskMean(fileList,toVecFunc,selVecFunc);
    C = diskCOV(fileList,U,toVecFunc,selVecFunc);
end