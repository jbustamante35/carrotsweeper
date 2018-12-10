function [area dB] = measureGerminate(Stack,BOX,mainBOX,numtoMeasure,disp)
    % get the centers of the kernels
    C = getCenters(Stack{1},mainBOX);
    %C = fliplr(C);
    clear BV
    %g1 = figure;
    %g2 = figure;
    % loop over the centers
    dB = {};
    for e = 1:size(C,1)
        dB{e} = {};
        try
           SKIP = 5;           
           A = [];
           for img = 1:SKIP:numtoMeasure
                tim = clock;
                I = getKernelImage(Stack{img},C(e,:),BOX);
                [MASK] = getObjectMask(I,10^5,1);
                dB{e}{end+1} = bwboundaries(MASK);
                A = [A sum(MASK(:))];
                dt = etime(clock,tim);
                fprintf(['Done with image in:' num2str(dt) '\n']);
           end
           area{e} = A;
           fprintf(['Done with:' num2str(e) '\n']);
        catch ME
            fprintf(['Error on:' num2str(e) '\n']);
        end
    end
    
end

%{
    
%}