classdef goT < matlab.mixin.Copyable
    properties
        %%%%%%%%%%
        % stack of position and direction
        position;
        direction;
        %%%%%%%%%%
        % parameters for nHood
        stepSize = .2;
        rho;
        rad;
        density;
        v1;
        v2;
        % nHood
        nHood_car;
        nHood_rad;
        W = [];
        %%%%%%%%%%
        % image to operate on
        image;
        %%%%%%%%%%
        % weighting function
        wFunc;
        
    end
    
    methods 
        function [obj] = goT()
            
        end
        
        %%%%%%%%%%
        % set position
        function [] = setPosition(obj,X)
            obj.position = X;
        end
        % set direction
        function [] = setDirection(obj,X)
            obj.direction = X;
        end
        
        
        %%%%%%%%%%
        % set rho
        function [] = setNhoodRho(obj,rho)
            obj.rho = rho;
        end
        % set rad
        function [] = setNhoodRad(obj,rad)
            obj.rad = rad;
        end
        % set nHood density
        function [] = setNhoodDensity(obj,density)
            obj.density = density;            
        end
        % generate the nHood
        function [] = generateH(obj)
            obj.v1 = linspace(0,obj.rho,obj.density(1));
            obj.v2 = linspace(-obj.rad,obj.rad,obj.density(2));
            [rho rad] = ndgrid(obj.v1,obj.v2);
            %obj.nHood_car = [rho(:)'.*cos(rad(:)');rho(:)'.*sin(rad(:)')];
            obj.nHood_car = [rho(:)'.*cos(rad(:)');rho(:)'.*sin(rad(:)')]';
            obj.nHood_rad = [rho(:)';rad(:)'];
            
        end
        % get nHood rot by direction
        function [rot] = getRot(obj,index)
            %rot = obj.direction(:,:,index)'*obj.nHood_car;
            rot = obj.nHood_car*obj.direction(:,:,index);
        end
        
        %%%%%%%%%%
        % set Image
        function [] = setImage(obj,I)
            obj.image = I;
        end
        
        %%%%%%%%%%
        % set weighting function
        function [] = setWfunction(obj,func)
            obj.wFunc = func;
        end
        
        %%%%%%%%%%
        % set Image
        function [] = setStepSize(obj,step)
            obj.stepSize = step;
        end
        
        
        
        % sample image at a point along the curve
        function [sam] = sampleImageAtPoint(obj,index)
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
        % sample the iamge over the  curve
        function [sam] = sampleImageAtCurve(obj,index)
            if nargin == 1;index(1) =1;index(2) = size(obj.position,2);end
            
            for e = index(1):index(2)
                sam(:,e) = obj.sampleImageAtPoint(e);
            end
        end
        % sample the curve at point
        function [sam] = sampleCurveAtPoint(obj,rho,index)
            tmpCurve = bsxfun(@minus,obj.position,obj.position(:,index));
            distanceToPoint = sum(tmpCurve.*tmpCurve,1).^.5;
            sidx = find(distanceToPoint < rho);
            sidx = sidx(sidx > index);
            sam = tmpCurve(:,sidx);
            sam = (obj.direction(:,:,index)'*sam);
            sam = goT.reparameterize(sam);
            sam = interp1(1:size(sam,2),sam',linspace(1,size(sam,2),rho))';
        end
        % sample the curve over the curve
        function [sam] = sampleCurveAtCurve(obj,rho)
            for e = 1:size(obj.position,2)-rho
                sam(:,:,e) = obj.sampleCurveAtPoint(rho,e);
            end
        end
        
        % arclength parameterize curve
        function [] = reparameterizeCurve(obj)
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
            nextPoint = obj.position(:,end) + obj.stepSize*vec;
            T = [vec]';
            N = [T(2) -T(1)];
            direction = ([T;N]);
        end
        
        function [] = walk(obj,steps,stepFunc)
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
                obj.position(:,end+1) = nP;
                obj.direction(:,:,end+1) = direc;
                %{
                % plot after step
                if nargin == 3;
                    obj.plotCurrentLocation(h);
                end
                %}
            end
        end
        
        function [] = walkUntil(obj,terminateFunc,stepFunc)
            while terminateFunc(obj.position);
                %{
                % plot before step
                if nargin == 3;
                    obj.plotCurrentLocation(h);
                end
                %}
                if nargin == 2
                    [nP direc] = obj.step();
                end
                obj.position(:,end+1) = nP;
                obj.direction(:,:,end+1) = direc;
                %{
                % plot after step
                if nargin == 3;
                    obj.plotCurrentLocation(h);
                end
                %}
            end
        end
        
        
        
        function [] = plotCurrentLocation(obj,h)
            plot(obj.position(1,end),obj.position(2,end),'b.');
            quiver(obj.position(1,end),obj.position(2,end),obj.direction(1,1,end),obj.direction(1,2,end),10,'Color','r');
            quiver(obj.position(1,end),obj.position(2,end),obj.direction(2,1,end),obj.direction(2,2,end),10,'Color','g');
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
            % construct bundle
            dT = gradient(curve);
            dL = sum(dT.*dT,1).^-.5;
            dT = bsxfun(@times,dL,dT);
            for e = 1:size(dT,2)
                TB(:,:,e) = [dT(:,e)';[dT(2,e) -dT(1,e)]];
            end
        end
        
        function [curve] = reparameterize(curve)
            dT = diff(curve,1,2);
            dT = sum(dT.*dT,1).^.5;
            L = cumsum([0 dT]);
            newL = linspace(0,L(end),round(L(end)));
            curve = interp1(L,curve',newL)';
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

%{
%%%%%%%%%%%%%%%%%%
% the width of the belief or the momentum is a function of percent 
% of cutoff
wsigma = .3;
maxWidth = 20*pi/180;
width = @(delta)maxWidth*normpdf(delta,0,wsigma)/normpdf(0,0,wsigma);
widthTest = linspace(0,1,100);
plot(widthTest,width(widthTest));
%%%%%%%%%%%%%%%%%%
% function for radial belief
k = 1.5;
alpha = 1;
scale = 10;
radial = @(x)wblpdf(x*scale^-1,alpha,k);
radialTest = linspace(0,20,100);
plot(radialTest,radial(radialTest));
%%%%%%%%%%%%%%%%%%
% function for angular belief
angle = @(x,delta)normpdf(x,0,width(delta));
%%%%%%%%%%%%%%%%%%
% belief function
func = @(rad,rho,delta)angle(rad,delta).*radial(rho);


pointDensity = [30 300];
RHO = 30;
RAD = pi;
[rho rad] = ndgrid(linspace(0,RHO,pointDensity(1)),linspace(RAD,-RAD,pointDensity(2)));
CAR = [rho(:)'.*cos(rad(:)');rho(:)'.*sin(rad(:)')];
PHI = [rho(:)';rad(:)'];
TEST = func(PHI(2,:),PHI(1,:),0);
TEST = reshape(TEST,pointDensity);
mesh(reshape(CAR(1,:),pointDensity),reshape(CAR(2,:),pointDensity),TEST);


I = double(imread(SET{10}{1}))/255;
[g1 g2] = gradient(I);
G = (g1.^2 + g2.^2).^.5;
[x y V] = impixel(I);
t1 = ba_interp2(g1,x,y);
t2 = ba_interp2(g2,x,y);
N = [t1 t2];
N = N / norm(N);
T = [-N(2) N(1)];
initD = [T;N];
T = goT();
T.setWfunction(func);
T.setNhoodRho(30);
T.setNhoodRad(pi);
T.setNhoodDensity(pointDensity);
T.generateH();
T.setPosition([x;y]);
T.setImage(G);
T.setDirection(initD);
T.walk(1000);
T.reparameterizeCurve();


PT = 100;
rho = 20;
for PT = 1:300
    samCurve = T.sampleCurveAtPoint(rho,PT);
    sam = T.sampleImageAtPoint(PT);
    sam = reshape(sam,[T.density]);
    rot = T.nHood_car;
    mesh(reshape(rot(1,:),[T.density]),reshape(rot(2,:),[T.density]),sam);
    hold on;
    plot(samCurve(1,:),samCurve(2,:),'k');
    view([0 -90]);
    drawnow
    hold off
end

%%
samCurve = T.sampleCurveAtCurve(20);
sam = T.sampleImageAtCurve();
sam = reshape(sam,[T.density size(sam,2)]);


figure
hold on
imshow(I,[]);hold on

T.plotCurrentPath();
T.plotFrameBundle();


%}