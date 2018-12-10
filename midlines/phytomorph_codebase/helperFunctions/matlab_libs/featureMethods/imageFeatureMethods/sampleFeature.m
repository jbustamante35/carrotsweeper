function [o] = sampleFeature(objStore,FileList,pointSet,para,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sample each point and obtain vector V
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %       objStore = Data base
    %       FileList = source of images
    %       pointSet = List of coordinates of corners
    %       para = Parameters for feature extractions and samplings. 
    %       disp = display. If disp = 1, the processed images are shown.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %       o = ?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        import phytoGraph_lite.Bpersist.Bos.newStore.*;
        import phytoGraph_lite.BBioData.phytoGeometry.tensors.*;
        of = 0;
            if (nargout == 1);
                of = 1;o.dVec = [];
                o.gVec = [];
                o.sz = [];
                o.pt=[];
                o.fn = {};
            end
        idisp = 1;
        manual = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create sample grid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % patch grid
        para = genDomains(para.dpara);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if disp;
            wbar = waitbar(0,'Sampling features...');
        end
        
        if idisp;
            dh = figure;
        end
        
        % loop images
        for i = 1:numel(pointSet)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read the image
            I = double(imread(FileList.get(i-1)));
            %subtract the background
            I = subtractBG(I);
               
            if idisp;
                    ptp = [];
                    vec1 = [];
                    vec2 = [];
                end
            featureSet = spfB(objStore);
            
            % loop points
                if manual;
                    pS = mySelector(I,pointSet{i});
                else;
                    pS = pointSet{i};
                end
            for j = 1:size(pS,1)
                curPoint = [pS(j,2) pS(j,1)];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Tin = [eye(2) curPoint'];
                [Tout Yout] = pbPCA((I),Tin,para{1}.D,para{2}.D);
                Yout = normalizeImage(Yout);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % featureBundle - for each point            
                nB = pointFeatureBundle(mat2bson(curPoint),mat2bson(Yout),mat2bson(Tout(:,1:2)),objStore);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % add bundle to feature set
                featureSet.insertfB(nB);
                if of;
                    o.dVec = [o.dVec Yout(:)];
                    o.gVec = [o.gVec;[i j]];
                    o.pt = [o.pt;curPoint];
                end
                
                if disp;
                    wbar = waitbar(j/numel(pointSet),
                    wbar,'Sampling features...');
                end
                
                if idisp;
                    ptp = [ptp;curPoint];
                    vec1 = [vec1;Tout(:,1)'];vec2 = [vec2;Tout(:,2)'];
                end
            end
            
            if of;
                o.fn{i} = FileList.get(i-1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mini_disp = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if idisp;
                figure(dh);
                imshow(I,[]);
                hold on;
                if ~isempty(ptp)
                    plot(ptp(:,1),ptp(:,2),'r.');                
                    quiver(ptp(:,1),ptp(:,2),vec1(:,1),vec1(:,2),.2,'g');
                    quiver(ptp(:,1),ptp(:,2),vec2(:,1),vec2(:,2),.2,'b');
                end
                EXP = maizeKN(char(FileList.get(i-1).getFullFileName()));
                title(num2str(EXP));
                hold off;
                drawnow;
            end
            if mini_disp
                % view the data                
                kineEditor(I,nB);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if of;
                o.sz = para{1}.sz;
            end
            
            if disp;
                fprintf(['Done with tangent space:' num2str(i) ':' num2str(numel(pointSet)) '\n']);
            end
            
            if idisp;
            end
        end
        if disp;
            close(wbar);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch ME
        Report = getReport(ME)
    end
end
