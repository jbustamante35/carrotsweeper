function [pointList] = obtainPointListFromPara(para)
    pointList(1,:) = para(2)*[cos(para(1)) sin(para(1))];
    pointList(2,:) = para(3)*[-sin(para(1)) cos(para(1))];
    pointList(3,:) = -para(2)*[cos(para(1)) sin(para(1))];
    pointList(4,:) = -para(3)*[-sin(para(1)) cos(para(1))];
end