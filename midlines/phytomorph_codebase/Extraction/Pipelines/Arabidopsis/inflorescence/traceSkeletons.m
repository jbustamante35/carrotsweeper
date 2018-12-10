function [PATHS] = traceSkeletons(fileList,idxL,startPoint,disp)
    I = imread(fileList{1});
    
    
    for s = 1:size(startPoint,1)
        for e = 1:numel(fileList)
            if disp
                tmpI = imread(fileList{e});
                imshow(tmpI,[]);
                hold on;
            end

            Z = zeros(size(I));
            fidx = find(idxL(:,2) == e);
            Z(idxL(fidx,1)) = 1;

            skel = [];
            SKEL = bwmorph(logical(Z),'skeleton',inf);
            [skel(:,1) skel(:,2)] = find(SKEL);
            midx = snapTo(skel,startPoint(s,:));
            skelStart = skel(midx,:);
            T = Radjacency(skel',2^.5+eps);


            EP = bwmorph(SKEL,'endpoints');
            ep = [];
            [ep(:,1) ep(:,2)] = find(EP);



            BP = bwmorph(SKEL,'branchpoints');
            bp = [];
            [bp(:,1) bp(:,2)] = find(BP);

            TP = [ep;bp];

            cnt = 1;
            for k = 1:size(TP,1)

                tidx = snapTo(skel,TP(k,:));
                [pidx,v] = dijkstra(T,midx, tidx);
                tpath = skel(pidx,:);
                if ~isinf(v) & v ~= 0
                    PATHS{s}{e}{cnt} = tpath;
                    cnt = cnt + 1;
                end

                %{ 
                imshow(tmpI,[]);
                hold on;
                plot(tpath(:,2),tpath(:,1),'m')
                hold off
                drawnow
                %}
            end
            %{
            pop = zeros(size(bp,1),1);
            for pth = 1:numel(PATHS{s}{e})
                for p = 1:size(bp,1)
                    if ~isempty(intersect(PATHS{s}{e}{pth},bp(p,:),'rows'))
                        pop(p) = pop(p) + 1;
                    end
                end
            end
            [pop,midx] = sort(pop,'descend');
            %[~,midx] = max(pop);
            sp = bp(midx(1:5),:);
            %}
            if disp
                plot(sp(:,2),sp(:,1),'r*');
                hold off
                drawnow
            end
        end
    end
end