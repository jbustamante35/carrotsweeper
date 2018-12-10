classdef arabidopsisImage < handle
    properties
        filename;
        imageData;
    end
    
    methods
        function [obj] = arabidopsisImage(filename)
            obj.filename = filename;
        end
        
        function [I] = imread(obj)
            obj.loadImage();
            I = obj.imageData;
        end
        
        function [] = loadImage(obj)
            if isempty(obj.imageData)
                obj.imageData = double(imread(obj.filename));
            end
        end
        
        function [curveSet] = extractClosedContours(obj,containerClass,contourLevel,filterLevel,filterLength)
            if nargin == 4;filterLength=[0 inf];end
            obj.loadImage();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reduce image to gray scale if color
            if size(obj.imageData,3) == 3
                obj.imageData = rgb2gray(obj.imageData/255);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % smooth the image
            tmp = imfilter(obj.imageData,fspecial('gaussian',5*filterLevel,filterLevel),'replicate');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obtain the gradient for contouring
            [g1 g2] = gradient(tmp);
            G = (g1.^2 + g2.^2).^.5;
            C = contourc(G,contourLevel);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % get the curve structure
            str = 1;
            c = 1;
            clear curve
            while str < size(C,2)
                ed = str + C(2,str);
                curve(c).level = C(1,str);
                curve(c).data = C(:,str+1:ed);
                curve(c).length = size(curve(c).data,2);
                c = c + 1;
                str = ed + 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % select closed curves and filter on length
            [curve] = selectClosedCurves(curve);
            ridx = [curve.length] < filterLength(1) & [curve.length] > filterLength(2);
            curve(ridx) = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % pour the curves into the given container class and assign the
            % file from which they came
            for e = 1:numel(curve)
                curveSet(e) = containerClass(curve(e).data);
                curveSet(e).filename = obj.filename;
            end
        end
        
        function [curveSet] = sampleGradient(obj,curves)
            obj.loadImage();
            [g1 g2] = gradient(obj.imageData);
            G = (g1.^2 + g2.^2).^.5;
            for e = 1:numel(curves)
                sample = ba_interp2(G,curves(e).segment(1,:),curves(e).segment(2,:));
                curves(e).totalGradient = mean(sample);
            end
        end
        
        function [I] = readCropped(obj,BOX)
            ROWS = round([BOX(1) BOX(1)+BOX(3)]);
            COLS = round([BOX(2) BOX(2)+BOX(4)]);
            imageData = double(imread(obj.filename,'PixelRegion',{COLS ROWS}));
            I = croppedMaizeImage( obj.filename);            
            I.cropBox = BOX;
            I.imageData = imageData;
        end
    end
end