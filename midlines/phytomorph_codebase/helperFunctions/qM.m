function [m] = qM(Observable,waveFunction1,waveFunction2)
    cnt = 1;
    for w2 = 1:size(waveFunction2,2)
        for w1 = 1:size(waveFunction1,2)
            waveFunction = waveFunction1(:,w1).*waveFunction2(:,w2);
            waveFunction = waveFunction / norm(waveFunction);
            for o = 1:numel(Observable)
                m(cnt,o) = waveFunction'*Observable{o}*waveFunction;
            end
            cnt = cnt +1;
        end
    end
end