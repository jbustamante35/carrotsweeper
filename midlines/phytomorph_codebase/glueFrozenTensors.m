function [T] = glueFrozenTensors(T1,T2)
    grade1 = T1(1,:);
    grade2 = T2(1,:);
    newGrade = grade1 + grade2;
    T = cat(1,newGrade,T1(2:end,:),T2(2:end,:));
end