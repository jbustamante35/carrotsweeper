function [D] = genD()
    vec = 2.^(0:8) + 1;
    for e = 1:9
        S1 = myGen3(vec(e));
        [q1 q2] = find(S1);
        for k = 1:9
            S2 = myGen3(vec(k));
            [w1 w2] = find(S2);
            d = ((q1 - w1)^2 + (q2 - w2)^2)^.5;
            if d <= 2^.5 & d ~=0
                D(e,k) = 1;
            else
                D(e,k) = 0;
            end
        end
    end
end