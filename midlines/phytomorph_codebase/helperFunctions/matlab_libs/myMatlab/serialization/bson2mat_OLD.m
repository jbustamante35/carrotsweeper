function [MAT] = bson2mat(BSON)
    tmp = myV(BSON);
    tmp.toMATLAB();
    MAT = tmp.d;
end