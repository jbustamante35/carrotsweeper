function [BSON] = mat2bson(MAT)
    tmp = myV(MAT);
    tmp.toBSON();
    BSON = tmp.d;
end