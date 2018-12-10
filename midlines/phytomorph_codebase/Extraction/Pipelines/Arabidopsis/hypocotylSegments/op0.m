function [] = op0(imageStack,oPath,disp)
    try
        numToProcess = numel(imageStack);
        numToProcess = min(numToProcess,50);
        %% remove not compliant images
        n = [];
        for e = 1:numToProcess
            try
                [p,tn,ex] = fileparts(imageStack{e});
                rm(e) = 0;
                if strcmp(tn(1),'.')
                    rm(e) = 1;
                end
            catch
                rm(e) = 1;
            end
        end
        imageStack(find(rm)) = [];

        %% sort the images
        n = [];
        for e = 1:numToProcess
            try
                [p,tn,ex] = fileparts(imageStack{e});
                n(e) = str2num(tn);    

            catch

            end
        end
        [~,sidx] = sort(n);
        imageStack = imageStack(sidx);

        %% read the images
        rm = [];
        for e = 1:numToProcess
            try
                I(:,:,e) = imread(imageStack{e});    
                rm(e) = 0;
            catch ME
                ME
                rm(e) = 1;
            end
        end
        I(:,:,find(rm)) = 1;

        %% binarize the images
        B = zeros(size(I,1),size(I,2),size(I,3));
        S = B;
        for e = 1:numToProcess
            tmp = double(I(:,:,e))/255;
            BK = imclose(tmp,strel('disk',61,0));
            BK = imfilter(BK,fspecial('gaussian',[51 51],9),'replicate');
            tmp = imfilter(tmp,fspecial('gaussian',[51 51],9),'replicate');
            tmp = tmp - BK;
            tmp = bindVec(tmp);
            level = graythresh(tmp);
            B(:,:,e) = double(tmp) < level;
            B(:,:,e) = bwareaopen(B(:,:,e),1000);
            B(:,:,e) = imfill(B(:,:,e),'holes');
            B(:,:,e) = imclearborder(B(:,:,e));
            S(:,:,e) = bwmorph(B(:,:,e),'skel',inf);
            tmp = S(:,:,e);
            for k = 1:20
                ep = bwmorph(tmp,'endpoints');
                tmp(find(ep)) = 0;
            end
            S(:,:,e) = tmp;
            e
        end


        %% get contours and midlines
        %dB = cell(numToProcess);
        P = cell(numToProcess,1);
        pidx = cell(numToProcess,1);
        loc = cell(numToProcess,1);

        for e = 1:numToProcess    
            dB{e} = bwboundaries(B(:,:,e));
            for c = 1:numel(dB{e})
                W = poly2mask(dB{e}{c}(:,2), dB{e}{c}(:,1), size(B,1),size(B,2));
                W = S(:,:,e).*W;
                [iy ix] = find(W);
                [out] = cwtK_closed_imfilter(dB{e}{c},{[51]});
                sig = [out.K;out.K;out.K];
                peaks = imdilate(-sig,strel('disk',50,0)) == -sig;
                str = size(out.K,1)+1;
                stp = str + size(out.K) -1;
                peaks = peaks(str:stp);
                pidx{e}{c} = find(peaks);
                values = out.K(pidx{e}{c});
                [~,sidx] = sort(values);
                pidx{e}{c} = pidx{e}{c}(sidx(1:2));
                loc{e}{c} = dB{e}{c}(pidx{e}{c},:);         
                DP = [fliplr(loc{e}{c})' [ix iy]'];
                T = Radjacency(DP,100);
                [path , pathcost]  = dijkstra(T , 1 , 2);
                P{e}{c} = DP(:,path);
                %fprintf(['Done with analysis: ' num2str(c) ' and image ' num2str(e) '\n']);
            end
        end


        if disp
            figure;
            %close all
            Z = zeros(size(B,1),size(B,2));
            for e = 1:numToProcess    
                RGB = cat(3,S(:,:,e),Z,B(:,:,e));        
                imshow(double(RGB),[]);
                hold on
                for c = 1:numel(dB{e})
                    plot(dB{e}{c}(:,2),dB{e}{c}(:,1),'g')
                    plot(loc{e}{c}(:,2),loc{e}{c}(:,1),'r*')
                    plot(P{e}{c}(1,:),P{e}{c}(2,:),'y')
                end
                hold off
                title(num2str(e))
                drawnow
            end
        end



        %% measure
        try
            L = [];
            WID = [];
            for img = 1:numToProcess
                tmp = bwdist(~B(:,:,img));
                for con = 1:numel(dB{img})

                    IDX = sub2ind(size(tmp),P{img}{con}(2,:),P{img}{con}(1,:));
                    values = tmp(IDX);
                    WID(img,con) = mean(values);
                    dL = diff(P{img}{con},1,2);
                    dL = sum(dL.*dL,1).^.5;
                    L(img,con) = sum(dL);
                    LEG{con} = num2str(con);
                end
                img
            end
        catch ME
                ME
        end


        %% plot measurements
        close all
        Z = zeros(size(B,1),size(B,2));
        plot(L)
        legend(LEG)
        figure
        plot(WID)
        legend(LEG)
        figure;
        e = 1;
        RGB = cat(3,S(:,:,e),Z,B(:,:,e));        
        %imshow(double(RGB),[]);
        I = imread(imageStack{1});
        %imshow(I,[]);
        % make display for condor
        image(I);
        hold on
        for c = 1:numel(dB{e})
            plot(dB{e}{c}(:,2),dB{e}{c}(:,1),'g')
            plot(loc{e}{c}(:,2),loc{e}{c}(:,1),'r*')
            plot(P{e}{c}(1,:),P{e}{c}(2,:),'y')
            uP = mean(dB{e}{c},1);
            text(uP(2),uP(1),num2str(c),'BackgroundColor',[1 1 1]);
        end
        hold off
        drawnow

        %% save the measurements
        [tp,tn,te] = fileparts(imageStack{1});

        fidx = strfind(tp,filesep);
        tnm = tp(fidx(end-2):end);
        tnm = strrep(tnm,filesep,'_');



        fileName = [oPath tnm '_frame1.tif'];    
        saveas(gca,fileName);

        %% plot measurements END
        close all
        Z = zeros(size(B,1),size(B,2));
        plot(L)
        legend(LEG)
        figure
        plot(WID)
        legend(LEG)
        figure;    
        I = imread(imageStack{end});
        imshow(I,[]);
        hold on
        for c = 1:numel(dB{end})
            plot(dB{end}{c}(:,2),dB{end}{c}(:,1),'g')
            plot(loc{end}{c}(:,2),loc{end}{c}(:,1),'r*')
            plot(P{end}{c}(1,:),P{end}{c}(2,:),'y')
            uP = mean(dB{end}{c},1);
            text(uP(2),uP(1),num2str(c),'BackgroundColor',[1 1 1]);
        end
        hold off
        drawnow

        %% save the measurements
        fileName = [oPath tnm '_frameEND.tif'];    
        saveas(gca,fileName);


        fileName = [oPath tnm 'length.csv'];
        csvwrite(fileName,L);

        fileName = [oPath tnm 'width.csv'];
        csvwrite(fileName,WID);
    catch ME
        close all
        [p,n,e] = fileparts(imageStack{1});
        p
    end
    close all
end