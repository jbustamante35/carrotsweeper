function [CI] = getPlantCells(FileList,GMModel,V,CL,oPath,disp)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the file parts, output Name (oName)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pth,nm,ext] = fileparts(FileList{1});
        sidx = strfind(pth,filesep);
        oName = pth((sidx(end-1)+1):end);
        oName = strrep(oName,filesep,'_');
        expand = 50;
        tform = [];
        resEstimate = NaN;

        
        I = imread(FileList{1});
        [nI,tform,resEstimate,tf] = getRectificedImage(I,GMModel,V,tform,resEstimate,disp);
        


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % average the stack
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cnt = 1;
        for e = 1:numel(FileList)
            % read image
            I = imread(FileList{e});
            if e == 1
                M = zeros(size(I));
            end
            % if day image
            if ~isNight(I,40)
                M = M + double(I)/255;
                cnt = cnt + 1;
            end
            e
        end
        M = M * cnt ^-1;
        L = labelImage(M*255,GMModel);
        [grandM,out] = label2Mask(M,L,V,CL,disp);
        masterBOX = mask2cellBox(grandM,10000,expand,disp,nI,CL);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each image
        R = {};
        for e = 1:numel(FileList)
            try
                % read image
                I = imread(FileList{e});
                % if day image
                if ~isNight(I,40)
                    close all force
                    % get rectified image
                    [nI,tform,resEstimate,tf] = getRectificedImage(I,GMModel,V,tform,resEstimate,disp);
                    if ~tf
                        % label image
                        L = labelImage(nI*255,GMModel);
                        % create masks from labeled image
                        [M,out] = label2Mask(nI,L,V,CL,disp);
                        % create crop boxes
                        %R{e} = mask2cellBox(M,10000,expand,disp,nI,CL);
                        % crop "cells"
                        CI{e} = cropCells(nI,masterBOX);
                        % crop labels
                        LI{e} = cropCells(L,masterBOX);
                        % get msg
                        msg{e} = getQRmsg(CI{e},LI{e},V);
                    end
                end

            catch ME

                R{e} = [];
                CI{e} = [];
                LI{e} = [];
                msg{e} = [];
               ME
            end

            [~,fn{e}] = fileparts(FileList{e});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        tmpI = imread(FileList{3});
        for e = 1:numel(CI)
            %{
            imshow(tmpI,[]);
            hold on
            %}
            nidx = [];
            if numel(CI{e}) == 9

                pointList = [];
                idList = (1:9)';

                for c = 1:numel(masterBOX)
                    pointList(c,:) = masterBOX(c).BoundingBox(1:2);
                end

                
                for row = 1:3
                    distance = sum(pointList.*pointList,2);
                    % find the upper left corner
                    [~,midx] = min(distance);
                    % find the first row
                    rowDistance = pointList(:,2) - pointList(midx,2);
                    [~,sidx] = sort(rowDistance);

                    [~,sidx1] = sort(pointList(sidx(1:3),1));


                    nidx = [nidx;idList(sidx(sidx1(1:3)))];


                    pointList(sidx(1:3),:) = [];
                    idList(sidx(1:3)) = [];

                end
                %R{e} = R{e}(nidx);
                R{e} = masterBOX(nidx);
                CI{e} = CI{e}(nidx);
                LI{e} = LI{e}(nidx);
                msg{e} = msg{e}(nidx);
                for pt = 1:numel(R{e})
                    %{
                    plot(R{e}(pt).BoundingBox(1),R{e}(pt).BoundingBox(2),'r*')
                    text(R{e}(pt).BoundingBox(1),R{e}(pt).BoundingBox(2),num2str(pt))
                    hold on
                    drawnow
                    %}
                end
                %hold off
                %waitforbuttonpress


            end
        end


        for e = 1:numel(LI)
            numCells(e) = numel(LI{e});
        end
        TnumCells = numCells;
        TnumCells(TnumCells==0) = [];
        N = mode(TnumCells);
        for tm = 1:numel(LI)
            for c = 1:numel(LI{tm})

                if ~isempty(LI{tm}) & numCells(tm) == N

                    if ~isempty(msg{tm})
                        MSG{c} = msg{tm}{c};
                    end

                    tmpM = label2Mask(LI{tm}{c},LI{tm}{c},V,CL,false);
                    plantMask = bwareaopen(tmpM(:,:,2),1000);
                    Area(tm,c) = sum(plantMask(:));
                end
            end
            tm
            numel(LI)
        end





        T = table;


        for t = 2:size(Area,1)
            T(t,'fileName') = {fn{t}};
        end
        %{
        for p = 1:size(Area,2)
            T(1,['plant_' num2str(p)]) = {MSG{p}};
        end
        %}

        for t = 2:size(Area,1)
            T(t,'fileName') = {fn{t}};
            for p = 1:size(Area,2)
                T(t,['plant_' num2str(p)]) = {num2str(Area(t,p))};
            end
        end


        writetable(T,[oPath oName '.csv']);

        %csvwrite([oPath oName '.csv'],Area);



        CMD = ['irsync -r -V /mnt/tetra/nate/overHeadReturn i:/iplant/home/cullenvens/overHeadReturn/']
        system(CMD);

        here = 1;
    catch ME
        getReport(ME)
        fprintf(['*************************\n'])
        fprintf(['*************************\n'])
         
    end
end