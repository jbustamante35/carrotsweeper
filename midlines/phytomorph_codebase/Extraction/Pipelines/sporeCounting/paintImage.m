function [out] = paintImage(I,patchStructure,oPath,fileName)
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    fprintf(['Start pixel painting\n']);tm = clock;
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    hyphaMask = zeros(size(I));
    sporeMask = zeros(size(I));

    for e = 1:numel(patchStructure)
        try
            hyphaMask(patchStructure(e).hyphaIDX) = 1;
            sporeMask(patchStructure(e).sporeIDX) = 1;
        catch
        end
    end
    out = flattenMaskOverlay(I,logical(hyphaMask),.4,'r');
    out = flattenMaskOverlay(out,logical(sporeMask),.3,'b');
    if ~isempty(oPath)
        imwrite(I,[oPath fileName '___original.jpg']);
        imwrite(out,[oPath fileName '___labled.jpg']);
    end
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    fprintf(['End pixel painting.' num2str(etime(clock,tm)) '\n']);
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
end