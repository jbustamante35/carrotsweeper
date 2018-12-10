function [out] = bson2myT(in)
    d = in.getData();
    d = bson2mat(d);
    sz = in.getSize();
end