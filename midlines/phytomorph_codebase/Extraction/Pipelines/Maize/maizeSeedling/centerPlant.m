function [tmpZ worked] = centerPlant(tmpI)
    worked = 1;
    tmpZ = zeros(1400,2100,'single');
    try

        tmpI((end-30):end,:) = [];

        %tmpI = logical(imresize(tmpI,.25));

        bot = tmpI(end-10:end,:);



        [dy dx] = find(bot);
        ux = round(mean(dx));

        [py px] = find(tmpI);

        px = px - ux;
        py = py - size(tmpI,1);

        px = px + 1050;
        py = py + 1400;

        idx = sub2ind(size(tmpZ),py,px);
        tmpZ(idx) = 1;
        tmpZ = imresize(tmpZ,[500 750]);
        tmpZ = single(tmpZ);
    catch
        worked = 0;
        tmpZ = imresize(tmpZ,[500 750]);
        tmpZ = single(tmpZ);
    end
end