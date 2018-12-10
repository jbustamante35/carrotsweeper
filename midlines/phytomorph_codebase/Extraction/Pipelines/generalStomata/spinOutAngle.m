function [para] = spinOutAngle(para)
    para(1) = atan2((para(2)/para(3))*sin(para(1)),cos(para(1)));
end