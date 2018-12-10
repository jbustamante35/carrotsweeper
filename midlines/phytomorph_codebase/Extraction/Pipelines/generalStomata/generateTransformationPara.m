function [PARA] = generateTransformationPara(E,O)
    % clicks gathered in counter-clockwise direction
    delta1 = E(1,:) - E(3,:);
    delta2 = E(2,:) - E(4,:);
    cp1 = E(3,:) + .5*delta1;
    cp2 = E(4,:) + .5*delta2;
    displacement = mean([cp1;cp2],1) - O;
    A1 = atan2(-delta1(2),-delta1(1));
    A2 = atan2(-delta2(2),-delta2(1));
    A = mean([A1;A2-pi/2]);
    PARA(1) = A;
    PARA(2) = norm(delta1)/2;
    PARA(3) = norm(delta2)/2;
    PARA(4) = displacement(1);
    PARA(5) = displacement(2);
end