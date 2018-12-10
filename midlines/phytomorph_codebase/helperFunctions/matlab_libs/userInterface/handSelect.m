function [data] = handSelect(FileList,data)
    
    model = 'null';
    data.codomain = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % loop over each image    
    for e = 1:FileList.size()
        
        
        selidx = (data.group(:,1) == e);
        
        pointSet = data.base(selidx,:);        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % guess
        [model grade] = initevalModel(data,model);
        subGrade = grade(selidx);
        [sGrade gidx] = sort(subGrade,'descend');
        gidx = gidx(1:data.target(e));
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % read in the image
        I = imread(FileList.get(e-1));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % get user input
        imshow(I,[]);hold on;
        plot(pointSet(:,1),pointSet(:,2),'r.');
        plot(pointSet(gidx,1),pointSet(gidx,2),'g*');
        [x1 x2 V] = impixel();
        selidx = zeros(size(pointSet,1),1);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % get user input
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % match to nearest point
        for f = 1:numel(x1)
            dist =  [];
            for g = 1:size(pointSet,1)
                dist(g) = norm(pointSet(g,:) - [x1(f) x2(f)]);
            end
            [j sidx] = min(dist);
            selidx(sidx) = 1;
            plot(pointSet(sidx,1),pointSet(sidx,2),'ko');
        end
        hold off;pause(.3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % create large label vector for codomain
        data.codomain = [data.codomain;selidx];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % build new model
        model = 'null';
        [model grade] = initevalModel(data,model);
        
        
        
        
        
        
    end
end