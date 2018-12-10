function [p] = toTM(numIn,numIn2)
    fidx = strfind(numIn2,':');
    numIn = num2str(numIn);
    t = datetime(str2num(numIn(1:4)),str2num(numIn(5:6)),str2num(numIn(7:8)),str2num(numIn2(1:(fidx(1)-1))),str2num(numIn2((fidx(1)+1):end)),0);
    p = posixtime(t);
end