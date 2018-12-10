classdef lManifold < handle
    
    properties
        rawX;       % raw X data        
        rawY;       % raw Y data
        simX;       % simX
        simY;       % simY
        kidx;       % groups
        nGroups;    % number of groups for learning
        initComp;   % init number of components pre kmeans clustering
        modelCompX; % modeling compX
        modelCompY; % modeling compY
        
        
        Ux;         % centers for X
        Uy;         % centers for Y
        Ex;         % basis for X
        Ey;         % basis for Y
        pVec;       % mapping from X -> Y
        corVecA;    % n corrlative vectors for the doamin
        corVecB;    % n corrlative vectors for the codoamin
        fit;        % information about fitting
    end
    
    
    methods
        function [obj] = lManifold()
            
        end
        %%%%%%%%%%%%%%%
        function [] = setX(obj,X)
            obj.rawX = X;
        end
        
        function [X] = getX(obj)
            X = obj.X;
        end
        
        function [] = clearX(obj)
            obj.setX([]);
        end
        
        function [] = addToX(obj,nX)
            obj.rawX = cat(1,obj.rawX,nX);
        end
        
        
        %%%%%%%%%%%%%%%
        function [] = setY(obj,Y)
            obj.rawY = Y;
        end
        
        function [Y] = getY(obj)
            Y = obj.Y;
        end
        
        function [] = clearY(obj)
            obj.setY([]);
        end
        
        function [] = addToY(obj,nY)
            obj.rawY = cat(1,obj.rawY,nY);
        end
        
        %%%%%%%%%%%%%%%        
        function [] = addXY(obj,nX,nY)
            addToX(obj,nX);
            addToY(obj,nY);
        end
        
        
        %%%%%%%%%%%%%%%
        function [] = setGroupN(obj,N)
            obj.nGroups = N;
        end
        
        function [] = setinitComp(obj,iC)
            obj.initComp = iC;
        end
        
        function [] = setmodelCompX(obj,mC)
            obj.modelCompX = mC;
        end
        
        function [] = setmodelCompY(obj,mC)
            obj.modelCompY = mC;
        end
        
        function [] = cleanData(obj)
            ridx = find(any(isnan(obj.rawY),2));
            obj.rawY(ridx,:) = [];
            obj.rawX(ridx,:) = [];
            
            ridx = find(any(isnan(obj.rawX),2));
            obj.rawY(ridx,:) = [];
            obj.rawX(ridx,:) = [];
        end
        
       
        function [] = learn(obj,N)            
            % decompose via linear
            [S C U E L ERR LAM] = PCA_FIT_FULL(obj.rawX,obj.initComp);
            % split into groups
            obj.kidx = kmeans(C,obj.nGroups);
            % loop over each group
            UQ = unique(obj.kidx);
            for u = 1:numel(UQ)                
                
                % set dataset X
                subX = obj.rawX(obj.kidx==UQ(u),:);
                [Sx Cx Ux Ex Lx ERRx LAMx] = PCA_FIT_FULL(subX,obj.modelCompX);
                obj.Ux(u,:) = Ux;
                obj.Ex(:,:,u) = Ex;
                %obj.Ex(:,:,u) = eye(size(Ex,1));
                Cx = PCA_REPROJ(subX,obj.Ex(:,:,u),obj.Ux(u,:));
                
                
                % sub dataset Y
                subY = obj.rawY(obj.kidx==UQ(u),:);                
                [Sy Cy Uy Ey Ly ERRy LAMy] = PCA_FIT_FULL(subY,obj.modelCompY);
                obj.Uy(u,:) = Uy;
                obj.Ey(:,:,u) = Ey;
                %obj.Ey(:,:,u) = eye(size(Ey,1));
                Cy = PCA_REPROJ(subY,obj.Ey(:,:,u),obj.Uy(u,:));
                
                
                
                % prediction
                %obj.pVec(u,:,:) = Cx'/Cy';
                [A,B,obj.fit.info{u},U,V,obj.fit.stats{u}] = canoncorr(Cx,Cy);
                obj.pVec(u,:,:) = A*inv(B);
                
                obj.corVecA{u} = A;
                obj.corVecB{u} = B;
            end
        end
        
        
        function [Y] = predict(obj,X)
            if isempty(X)
                Y = X;
                return
            end
            % for each vector in X
            for e = 1:size(X,1)
                % find the center each dataset maps to
                deltaCenter = bsxfun(@minus,obj.Ux,X(e,:));
                deltaCenter = sum(deltaCenter.*deltaCenter,2);
                [J,midx] = min(deltaCenter);
                % project data into the frame at that center
                %newX = X(e,:) - obj.Ux(midx,:);
                %newX = newX*obj.Ex(:,:,midx);
                newX = PCA_REPROJ(X(e,:),obj.Ex(:,:,midx),obj.Ux(midx,:));
                % create prediction
                predictY = newX*shiftdim(obj.pVec(midx,:,:),1);
                % back project and offset the mean
                %predictY = obj.Ey(:,:,midx)*predictY';
                %predictY = predictY + obj.Uy(midx,:)';
                predictY = PCA_BKPROJ(predictY,obj.Ey(:,:,midx),obj.Uy(midx,:));
                % store the result
                Y(e,:) = predictY';
            end
        end
        
    end
end