function [FUN GRID IG IGM] = ringSample(imgToSample,M,dB,R,NI,disp)
    try
        % if display
        if disp
            imshow(imgToSample,[])
            hold on
        end
        % init the vars out
        FUN = [];
        GRID = [];
        IG = {};
        % for each contour
        for e = 1:numel(dB)
            fprintf(['start:' num2str(e) ';' num2str(numel(dB)) '\n']);
            % crop the mask and the image
            IG{e} = imcrop(imgToSample,round(R(e).BoundingBox));
            Z = logical(imcrop(M,round(R(e).BoundingBox)));
            Z = bwlarge(Z);
            IGM{e} = Z;
            % take steps
            for d = 1:60
                % get the contour
                d2B = bwboundaries(logical(Z));
                % if disp
                if disp                
                    plot(d2B{1}(:,2),d2B{1}(:,1),'r')
                    drawnow
                end
                % sample along the contour
                tmp = [];
                for p = 1:size(d2B{1},1)
                    for k = 1:3
                        tmp(p,k) = IG{e}(d2B{1}(p,1),d2B{1}(p,2),k);
                    end                
                end
                % interp the same number of points along the contour and grid
                GRID(:,:,d,e) = interp1(1:size(d2B{1}),d2B{1},linspace(1,size(tmp,1),NI));
                FUN(:,:,d,e) = interp1(1:size(tmp,1),tmp,linspace(1,size(tmp,1),NI));
                % step in
                Z = imerode(logical(Z),strel('disk',2,0));
                Z = bwlarge(Z);
            end
            fprintf(['end:' num2str(e) ';' num2str(numel(dB)) '\n']);
        end
    catch ME
        ME
        FUN = [];
        GRID = [];
        IG = {};
    end
end