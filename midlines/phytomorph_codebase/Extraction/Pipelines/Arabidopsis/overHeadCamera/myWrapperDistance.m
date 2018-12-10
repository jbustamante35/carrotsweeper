function [dis] = myWrapperDistance(data,loc,max)
    mdata2 = data+max;
    mdata1 = data-max;
    dd = [(mdata1-loc) (data-loc) (mdata2-loc)];
    sd = sign(dd);
    [dis,idx] = min(abs(dd),[],2);
    idx = sub2ind(size(sd),(1:size(sd,1))',idx);
    dis = sd(idx).*dis;
end