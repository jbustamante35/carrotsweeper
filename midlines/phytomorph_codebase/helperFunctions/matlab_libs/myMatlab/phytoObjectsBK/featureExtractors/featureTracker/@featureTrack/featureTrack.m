classdef featureTrack  < phytoFuncSequence
    % isa sequence of calls 1:featureMap_atP -->
    properties
    
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = featureTrack(varargin)
            obj.para = [];
            obj.notes = [];
            obj.func = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set width and height of crop box
        function [] = setWH(obj,w,h)
            obj.para.w = w;
            obj.para.h = w;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % attach the feature function
        function [] = attachFeatureFunction(obj,f)
            obj.putFunc(f,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % attach the point selector
        function [] = attachPointSelector(obj,f)
            obj.putFunc(f,2);            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % track
        function [] = track(obj,X,P)
            try
                nT = P.nTrials();
                for file = 1:(numel(X)-1)
                    for pt = 1:nT
                        
                        % current point
                        curPoint = P(:,file,pt); 
                        
                        
                        % current cropBox
                        cropBox = phytoFunc.constructCropBox_atP(curPoint,obj.para.w,obj.para.h);
                        % compute feature map @ P
                        fM = obj.func{1}.computeFeatureMap_atP(X{file+1},cropBox);                    
                        % extract points
                        pL = obj.func{2}.extractPoints(fM);                    
                        
                        % extract matrix
                        pL = pL.d;
                        
                        % get center point of patch
                        sz = size(fM);
                        centerPoint = (fliplr((sz(1:2)-1)/2) + 1)';
                        % center points
                        pL = bsxfun(@minus,pL,centerPoint);
                        % delta vector
                        deltaV = obj.distanceMeasure([0;0],pL);
                        
                        % set next point
                        P(:,file+1,pt) = curPoint + deltaV;
                        fprintf(['point@' num2str(pt) ':' num2str(nT) ':' 'frame@' num2str(file) ':' num2str(numel(X)) '\n']);
                    end                
                end
                % set view
                para.type = 'phytoCurve';
                para.color = 'r';                    
                P.setView(para);
            catch ME
                ME
            end
        end
        
       
        
    end
    
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % distance measure
        function [p] = distanceMeasure(x,Y)
            if ~isempty(Y)
                for e = 1:size(Y,2)
                    delta = x - Y(:,e);
                    dist(e) = norm(delta);
                end
                [J idx] = min(dist);
                p = Y(:,idx);
            else
                p = x;
            end
        end
    end
    
    methods
        %%%%%%%%%%%%
        % eval chain
        function res = e(obj,X)
            
            obj.f
        end
        
    end
end