function [D] = modelGeneration(dataFiles,h)
    D = [];
    for i = 1:numel(dataFiles)
        % load
        dataS = genLoader(dataFiles{i},i);
        S = dataS.imageStack;
        I = imread(S{1}.fileName);
        J = zeros([size(I) S.size()]);
        G = zeros(size(J));
        parfor j = 1:61
            I = imread(S{j}.fileName);
            J(:,:,j) = I;
            [g1 g2] = gradient(double(I));
            g = (g1.^2 + g2.^2).^.5;
            G(:,:,j) = g;
        end

        sG = std(G,1,3);
        sJ = std(double(J),1,3);

        sJ = sJ.*sG;

        sJ = normalize(sJ);
        thresh = graythresh(sJ);
        MASK = sJ > thresh;
        MASK = imclose(MASK,strel('disk',21));
        MASK = imerode(MASK,ones(5));

        MASK = bwareaopen(MASK,100);
        CC = bwconncomp(MASK,8);
        R = regionprops(MASK,'Centroid');



        clickPos = [R.Centroid];
        clickPos = reshape(clickPos,[2 numel(R)])';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % snap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sidx = [];
        for k = 1:size(clickPos,1)
            dist = [];
            for l = 1:size(dataS.gamma,1)
                delta = clickPos(k,:) - dataS.gamma(l,:);
                dist(l) = norm(delta);
            end
            [J sidx(k)] = min(dist);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for view points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % draw data for user
        close all
        I = myReader(dataS.imageStack{1}.fileName,'toGray',1);                
        imshow(I);
        hold on;
        nProps.Color = 'g';
        dataS.gD.pointList{1}.view(h,nProps);   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view graphs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for u = 1:numel(sidx)
            vProps.Color = 'r';
            dataS.G.sG{sidx(u)}.N{1}{1}.view(h,vProps);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        drawnow
        D = [D;dataS.V(sidx,:)];
    end
end