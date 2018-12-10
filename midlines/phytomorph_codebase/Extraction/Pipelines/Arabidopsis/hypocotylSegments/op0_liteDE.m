function [] = op0_liteDE(imageStack,oPath,disp)
    try
        [p,n,e] = fileparts(imageStack{1});
        mkdir(oPath);
        % make the num to process small for testing
        numToProcess = numel(imageStack);
        %numToProcess = min(numToProcess,10);
        % remove and sort imageStack
        [imageStack] = removeANDsort(imageStack);

        for e = 1:numToProcess
            [L{e} WID{e} dB{e} P{e} loc{e}] = opOnSingleImage(imageStack{e},1);
        end
        for e = 1:numel(L)
            for c = 1:numel(L{e})
                wL(e,c) = L{e}(c);
                wW(e,c) = WID{e}(c);
            end
        end

        

        %% save the measurements
        [tp,tn,te] = fileparts(imageStack{1});

        %fidx = strfind(tp,filesep);
        %tnm = tp(fidx(end-2):end);
        %tnm = strrep(tnm,filesep,'_');
        figure;    
        I = imread(imageStack{1});
        I = cat(3,I,I,I);
        image(I);
        %imshow(I,[]);
        hold on
        for c = 1:numel(dB{1})
            plot(dB{1}{c}(:,2),dB{1}{c}(:,1),'g')
            plot(loc{1}{c}(:,2),loc{1}{c}(:,1),'r*')
            plot(P{1}{c}(1,:),P{1}{c}(2,:),'y')
            uP = mean(dB{1}{c},1);
            text(uP(2),uP(1),num2str(c),'BackgroundColor',[1 1 1]);
        end
        hold off
        drawnow

        fileName = [oPath p '_' tn '_frame1.tif'];    
        saveas(gca,fileName);
        
        %% plot measurements END
        close all
        %Z = zeros(size(B,1),size(B,2));
        %{
        plot(L)
        legend(LEG)

        figure
        plot(WID)
        legend(LEG)
        %}

        figure;    
        I = imread(imageStack{end});
        I = cat(3,I,I,I);
        image(I);
        %imshow(I,[]);
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
        fileName = [oPath p '_' tn '_frameEND.tif'];    
        saveas(gca,fileName);


        fileName = [oPath p '_' tn '_length.csv'];
        csvwrite(fileName,wL);

        fileName = [oPath p '_' tn '_width.csv'];
        csvwrite(fileName,wW);
        close all
    catch ME
        close all
        getReport(ME)
        [p,n,e] = fileparts(imageStack{1});
        p
    end
    close all
end

%{
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
            % read, double, normalize
            tmp = double(I(:,:,e))/255;
            % close the image
            BK = imclose(tmp,strel('disk',61,0));
            % filter and smooth the closed image
            BK = imfilter(BK,fspecial('gaussian',[51 51],9),'replicate');
            % filter the orginal
            tmp = imfilter(tmp,fspecial('gaussian',[51 51],9),'replicate');
            % subtract off the background
            tmp = tmp - BK;
            % normalize
            tmp = bindVec(tmp);
            % threshold value
            level = graythresh(tmp);
            % threshold image
            B(:,:,e) = double(tmp) < level;
            % remove greater than 1000 objects
            B(:,:,e) = bwareaopen(B(:,:,e),1000);
            % fill holes
            B(:,:,e) = imfill(B(:,:,e),'holes');
            % clear border
            B(:,:,e) = imclearborder(B(:,:,e));
            % get skeleton
            S(:,:,e) = bwmorph(B(:,:,e),'skel',inf);
            % store skeleton in tmp
            tmp = S(:,:,e);
            % process skeleton 
            for k = 1:20
                ep = bwmorph(tmp,'endpoints');
                tmp(find(ep)) = 0;
            end
            % store skeleton
            S(:,:,e) = tmp;
            e
        end


        %% get contours and midlines
        %dB = cell(numToProcess);
        P = cell(numToProcess,1);
        pidx = cell(numToProcess,1);
        loc = cell(numToProcess,1);

        for e = 1:numToProcess  
            % get the boundary of the mask
            dB{e} = bwboundaries(B(:,:,e));
            for c = 1:numel(dB{e})
                % get the poly2 mask
                W = poly2mask(dB{e}{c}(:,2), dB{e}{c}(:,1), size(B,1),size(B,2));
                % select the skeleton based on mask
                W = S(:,:,e).*W;
                % find the skeleton
                [iy ix] = find(W);
                % measure the curvature
                [out] = cwtK_closed_imfilter(dB{e}{c},{[51]});
                % stack the Kurvature
                sig = [out.K;out.K;out.K];
                % dilate the stacked curvature
                peaks = imdilate(-sig,strel('disk',50,0)) == -sig;
                % get the str
                str = size(out.K,1)+1;
                % get the stop
                stp = str + size(out.K) -1;
                % clip out the potential peaks
                peaks = peaks(str:stp);
                % find the peaks
                pidx{e}{c} = find(peaks);
                % get the kurvature values 
                values = out.K(pidx{e}{c});
                % sort the peaks
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
%}