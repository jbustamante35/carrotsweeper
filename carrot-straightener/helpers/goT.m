%% goT: class to hold data for tracking curve data
% change

classdef goT < matlab.mixin.Copyable
    properties
        % stack of position and direction
        position
        direction
        % parameters for nHood
        STEP_SIZE = 0.2
        rho
        rad
        density
        v1
        v2
        
        % nHood
        nHood_car
        nHood_rad
        W = []
        
        % image to operate on
        image
        
        % weighting function
        wFunc
    end
    
    methods
        function [obj] = goT()
            
        end
        
        %% set position
        function [] = setPosition(obj,X)
            obj.position = X;
        end
        
        %% set direction
        function [] = setDirection(obj,X)
            obj.direction = X;
        end
        
        function [] = setNhoodRho(obj,rho)
            %% set rho
            obj.rho = rho;
        end
        
        function [] = setNhoodRad(obj,rad)
            %% set rad
            obj.rad = rad;
        end
        
        function [] = setNhoodDensity(obj,density)
            %% set nHood density
            obj.density = density;
        end
        
        function [] = generateH(obj)
            %% generate the nHood
            obj.v1 = linspace(0,obj.rho,obj.density(1));
            obj.v2 = linspace(-obj.rad,obj.rad,obj.density(2));
            [rho rad] = ndgrid(obj.v1,obj.v2);
            %obj.nHood_car = [rho(:)'.*cos(rad(:)');rho(:)'.*sin(rad(:)')];
            obj.nHood_car = [rho(:)'.*cos(rad(:)');rho(:)'.*sin(rad(:)')]';
            obj.nHood_rad = [rho(:)';rad(:)'];
            
        end
        
        function [rot] = getRot(obj,index)
            %% get nHood rot by direction
            %rot = obj.direction(:,:,index)'*obj.nHood_car;
            rot = obj.nHood_car*obj.direction(:,:,index);
        end
        
        function [] = setImage(obj,I)
            %% set Image
            obj.image = I;
        end
        
        function [] = setWfunction(obj,func)
            %% set weighting function
            obj.wFunc = func;
        end
        
        function [] = setStepSize(obj,step)
            %% set Image
            obj.STEP_SIZE = step;
        end
        
        function [sam] = sampleImageAtPoint(obj,index)
            %% sample image at a point along the curve
            if nargin < 2;index = size(obj.position,2);end
            % rotate to last frame
            %tmpH = obj.direction(:,:,index)'*obj.nHood_car;
            tmpH = obj.nHood_car*obj.direction(:,:,index);
            % displace to last position
            %tmpH = bsxfun(@plus,tmpH,obj.position(:,index));
            tmpH = bsxfun(@plus,tmpH,obj.position(:,index)');
            % interpolate
            %sam = ba_interp2(obj.image,tmpH(1,:),tmpH(2,:));
            sam = ba_interp2(obj.image,tmpH(:,1),tmpH(:,2));
        end
        
        function [sam] = sampleImageAtCurve(obj,index)
            %% sample the iamge over the  curve
            if nargin == 1;index(1) =1;index(2) = size(obj.position,2);end
            
            for e = index(1):index(2)
                sam(:,e) = obj.sampleImageAtPoint(e);
            end
        end
        
        function [sam] = sampleCurveAtPoint(obj,rho,index)
            %% sample the curve at point
            tmpCurve = bsxfun(@minus,obj.position,obj.position(:,index));
            distanceToPoint = sum(tmpCurve.*tmpCurve,1).^.5;
            sidx = find(distanceToPoint < rho);
            sidx = sidx(sidx > index);
            sam = tmpCurve(:,sidx);
            sam = (obj.direction(:,:,index)'*sam);
            sam = goT.reparameterize(sam);
            sam = interp1(1:size(sam,2),sam',linspace(1,size(sam,2),rho))';
        end
        
        function [sam] = sampleCurveAtCurve(obj,rho)
            %% sample the curve over the curve
            for e = 1:size(obj.position,2)-rho
                sam(:,:,e) = obj.sampleCurveAtPoint(rho,e);
            end
        end
        
        function [] = reparameterizeCurve(obj)
            %% arclength parameterize curve
            dT = diff(obj.position,1,2);
            dT = sum(dT.*dT,1).^.5;
            L = cumsum([0 dT]);
            newL = linspace(0,L(end),round(L(end)));
            % for position
            sig = obj.position';
            obj.position = interp1(L,sig,newL)';
            % for frame bundle
            sz = size(obj.direction);
            sig = reshape(obj.direction,[sz(1)*sz(2) sz(3)]);
            sig = interp1(L,sig',newL);
            obj.direction = reshape(sig',[sz(1) sz(2) numel(newL)]);
        end
        
        
        function [nextPoint direction] = step(obj)
            %% step function
            if isempty(obj.W)
                % eval weight function
                obj.W = obj.wFunc(obj.nHood_rad(2,:),obj.nHood_rad(1,:),0)';
            end
            
            % sample the image
            samp = obj.sampleImageAtPoint();
            
            % multiply
            f = obj.W.*samp;
            
            % reshape
            f = reshape(f,obj.density);
            
            % make choice
            H = sum(f,1);
            [J,sidx] = max(H);
            theta = obj.v2(sidx);
            vec = [cos(theta);sin(theta)];
            vec = obj.direction(:,:,end)*vec;
            
            % get next point and direction
            nextPoint = obj.position(:,end) + obj.STEP_SIZE * vec;
            T = [vec]';
            N = [T(2) -T(1)];
            direction = ([T;N]);
        end
        
        function [] = walk(obj, steps, stepFunc)
            %% walk function
            for e = 1:steps
                %{
                % plot before step
                if nargin == 3;
                    obj.plotCurrentLocation(h);
                end
                %}
                
                if nargin == 2
                    [nP direc] = obj.step();
                end
                
                obj.position(:,end+1)    = nP;
                obj.direction(:,:,end+1) = direc;
                %{
                % plot after step
                if nargin == 3;
                    obj.plotCurrentLocation(h);
                end
                %}
            end
        end
        
        function [] = walkUntil(obj, terminateFunc, stepFunc)
            %% walk function until stop point
            while terminateFunc(obj.position)
                
                % plot before step
                %                 if nargin == 3
                %                     obj.plotCurrentLocation(h);
                %                 end
                
                if nargin == 2
                    [nP , direc] = obj.step();
                end
                
                obj.position(:, end+1)     = nP;
                obj.direction(:, :, end+1) = direc;
                
                % plot after step
                %                 if nargin == 3
                %                     obj.plotCurrentLocation(h);
                %                 end
                
            end
        end
        
        
        
        function [] = plotCurrentLocation(obj, h)
            plot(obj.position(1,end), obj.position(2,end), 'b.');
            
            quiver(obj.position(1,end), obj.position(2,end), ...
                obj.direction(1,1,end), obj.direction(1,2,end), ...
                10, 'Color', 'r');
            
            quiver(obj.position(1,end), obj.position(2,end), ...
                obj.direction(2,1,end), obj.direction(2,2,end), ...
                10, 'Color', 'g');
        end
        
        function [] = plotCurrentPath(obj)
            plot(obj.position(1,:),obj.position(2,:),'r');
        end
        
        function [] = plotFrameBundle(obj)
            quiver(obj.position(1,:),obj.position(2,:),squeeze(obj.direction(1,1,:))',squeeze(obj.direction(1,2,:))',10,'Color','r');
            quiver(obj.position(1,:),obj.position(2,:),squeeze(obj.direction(2,1,:))',squeeze(obj.direction(2,2,:))',10,'Color','g');
        end
        
    end
    
    methods (Static)
        function [TB] = generateTangentBundle(curve)
            %% construct bundle
            dT = gradient(curve);
            dL = sum(dT .* dT, 1).^-0.5;
            dT = bsxfun(@times, dL, dT);
            
            for e = 1:size(dT,2)
                TB(:,:,e) = [dT(:,e)' ; [dT(2,e) -dT(1,e)]];
            end
        end
        
        function [curve] = reparameterize(curve)
            dT    = diff(curve, 1, 2);
            dT    = sum(dT .* dT, 1).^0.5;
            L     = cumsum([0 dT]);
            newL  = linspace(0, L(end), round(L(end)));
            curve = interp1(L, curve', newL)';
        end
        
        function [] = viewData(curveData,imageData)
            for i = 1:size(curve,1)
                
                mesh(reshape(rot(1,:),[T.density]),reshape(rot(2,:),[T.density]),sam);
                hold on;
                plot(samCurve(1,:),samCurve(2,:),'k');
            end
            
        end
    end
end
