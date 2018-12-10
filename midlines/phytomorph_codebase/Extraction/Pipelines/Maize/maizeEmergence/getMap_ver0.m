function [sIDX] = getMap_ver0(tMASK,rec,OFFSET)
        try


        disp = 0;
        R = regionprops(tMASK,'Centroid');
        for e = 1:numel(R)
            delta(e) = sum(R(e).Centroid.^2);
        end
        %imshow(tMASK,[]);
        cnt = 1;
        MN = 2;
        [~,Sidx] = sort(delta);
        [~,sidx] = min([R(Sidx(1)).Centroid(2);R(Sidx(2)).Centroid(2)]);
        Sidx = Sidx(sidx);
        CM = [];
        CM = [CM;R(Sidx(end)).Centroid];
        R(Sidx) = [];
        while ~isempty(R)
            if disp
                imshow(tMASK,[]);
                hold on
                plot(CM(end,1),CM(end,2),'r*')
                drawnow
            end

            % check for those to the right
            delta = [];
            sel = [];
            dis = [];
            for e = 1:numel(R)
                dis(e,:) = R(e).Centroid - CM(end,:);
                delta(e) = sum(dis(e,:).*dis(e,:),2);
                if dis(e,1) > 0
                    sel(e) = 1;
                else
                    sel(e) = 0;
                end
            end
            fidx = find(sel);


            if (rem(size(CM,1),12) ~= 0 & cnt ~= 7)
                [J,sidx] = sort(delta(fidx));
                try
                next = fidx(sidx(1:min(MN,numel(delta))));
                catch ME
                    ME
                end
                if disp
                    for d = 1:numel(next)
                        plot(R(next(d)).Centroid(1),R(next(d)).Centroid(2),'go')
                    end
                end
                
                [J,nn] = sort(dis(next,2));
                if disp
                    plot(R(next(nn(1))).Centroid(1),R(next(nn(1))).Centroid(2),'b*');
                    drawnow
                    %waitforbuttonpress
                end

                next = next(nn(1));
            else
                
                if cnt == 13
                    MN = 1;
                end
                delta = [];
                
                if rem(size(CM,1),12) ~= 0
                    tmpCP = CM(end,:);
                else
                    tmpCP = [0 0];
                end
                
                for e = 1:numel(R)
                    delta(e) = sum((R(e).Centroid - tmpCP).^2);
                end
                cnt;
                [~,next] = min(delta);
                next = next(1);
                close all
                if rem(size(CM,1),12) == 0
                    cnt = cnt + 1;
                end
            end
        %hold off

            Sidx = [Sidx;next];
            CM = [CM;R(Sidx(end)).Centroid];



            R(next) = [];
            %drawnow
        end

        sIDX = CM;


        %{
        iVal = sum(tMASK,2) + rand(size(tMASK,1),1);
        iVal = imfilter(iVal,ones(75,1));
        rowIDX = find(iVal == imdilate(iVal,ones(150,1)));
        [J,sidx] = sort(iVal(rowIDX),'descend');
        rowIDX = rowIDX(sidx(1:14));
        rowIDX = sort(rowIDX);
        %}

        %{
        plot(iVal,'b');
        z = zeros(size(iVal));
        z(rowIDX) = 1;
        hold on
        plot(z*max(iVal),'r');
        %}
        %{
        sIDX = [];
        for e = 1:numel(rowIDX)
            z = zeros(size(tMASK,1),1);
            z(rowIDX(e)) = 1;
            z = imdilate(z,ones([200 1]));
            z = repmat(z,[1 size(tMASK,2)]);
            for i = 1:numel(R)
                if z(R(i).Centroid(2),R(i).Centroid(1))
                    sIDX = [sIDX;[R(i).Centroid]];
                end
            end
        end
        %}
        %% make labels
        LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
        NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
        HL = {'d' 'p'};
        LABELS = {};
        for e1 = 1:numel(HL)
            for e2 = 1:numel(LL)
                for e3 = 1:numel(NL)
                    LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
                end
            end
        end
        
        [sIDX(:,1) sIDX(:,2)] = tforminv(rec,sIDX(:,1)+OFFSET(1),sIDX(:,2)+OFFSET(2));
        MAP.centers = sIDX;
        MAP.labels = LABELS;
        
        %{
        close all
        imshow(tMASK,[]);
        hold on
        plot(sIDX(:,1),sIDX(:,2))
        %}
        catch ME
            getReport(ME)
            here = 1;
        end
end