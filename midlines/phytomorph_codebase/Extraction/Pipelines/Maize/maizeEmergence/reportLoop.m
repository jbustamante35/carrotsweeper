function [] = reportLoop(message,loopVar,modN)
    if loopVar == 1
            fprintf(['\n' message '\n']);
    end
    fprintf(['.']);
    if mod(loopVar,modN) == 0;
        fprintf(['\n'])
    end
end