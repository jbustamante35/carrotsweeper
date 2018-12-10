function [dE,C] = recursiveTC(data,U1,E1,U2,E2,modelStruct)

    [dE1,C1,errDis1] = applyTC(data,U1{1},E1{1},U2{1},E2{1},modelStruct{1});
    errDis1 = reshape(errDis1,size(data));
    [dE2,C2,errDis2] = applyTC(errDis1,U1{2},E1{2},U2{2},E2{2},modelStruct{2});
    C = [C1;C2];
    dE = [dE1;dE2];
    errDis = 0;
end