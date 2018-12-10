function[] = modelSweeper(M,label,oPath)
    numToSweep = numel(M.V);
    for c = 1:numToSweep
        figure;
        tmpU = zeros(numel(M.V));
        L = linspace(-M.V(c).^.5,M.V(c).^.5,7);
        for l = 1:numel(L)
            tmpC = tmpU;
            tmpC(c) = L(l);
            sweepC = PCA_BKPROJ_T(tmpC,M.E,M.U);
            plot(sweepC);
            hold all
        end
        title([label 'PC-sweep' num2str(c)]);
        saveas(gca,[oPath label 'PC-sweep' num2str(c) '.tif']);
        %waitforbuttonpress
        %close all
    end
end
