function [] = myClickTracker(varargin)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define inpath
    if nargin == 0
        inFilePath = uigetdir();
        inFilePath = [inFilePath filesep];
    else
        inFilePath = varargin{1};
    end
    % define outpath
    %oPath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction_14.03.05/';
    %mkdir(oPath)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scan for new images
    FileList = {};
    FileExt = {'tiff','TIF','tif'};
    verbose = 1;
    SET = sdig(inFilePath,FileList,FileExt,1); 
    
    
    cornerViewMovie = 0;
    if cornerViewMovie
        mov = figure;
        for e = 1:numel(SET)
            for f = 1:numel(SET{e})
                I = imread(SET{e}{f});
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find all corner(s) in the first image
                cornersFeatures = cornerFeatureMap(5);
                cornerFeatures = cornersFeatures.computeFeatureMap(double(I));
                allCorners = simplePointExtractor(20,.08);
                allCorners = allCorners.extractPoints(cornerFeatures);
                
                
                imshow(I,[])
                hold on
                allCorners.view(gca,[]);
                drawnow
                hold off
            end
            
        end
    end
    
    
    for e = 1:numel(SET)
        I = imread(SET{e}{1});
        [pth nm ext] = fileparts(SET{e}{1});
        
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find all corner(s) in the first image
        cornersFeatures = cornerFeatureMap(5);
        cornerFeatures = cornersFeatures.computeFeatureMap(double(I));
        allCorners = simplePointExtractor(20,.08);
        allCorners = allCorners.extractPoints(cornerFeatures);
        
        
        [r c V] = impixel(I);
        
        tracePath = myHS_X('phytoAcurve');
        for g = 1:numel(r)
            dist = @(x)sum(([r(g) c(g) 1]-x(1)).^2,2).^.5;
            dist = allCorners.distrib(dist);
            [JUNK midx] = min(cell2mat(dist(:)));
            tracePath{g} = phytoAcurve.contructFromApoint(allCorners{midx});
        end
        myQuickTracker(SET{e},tracePath);
        imshow(I,[])
        hold on
        tracePath.view(gca,[]);
        waitforbuttonpress
        L = tracePath.iLength();
        length = [];
        for i = 1:numel(L)
            length = [length L{i}(1:250)];
        end
        close all        
        figure;
        plot(length);
        waitforbuttonpress
        csvwrite([pth 'results.csv'],length);
    end
end
%{
    myClickTracker('/mnt/spaldingdata/nate/Nate/');
    myClickTracker('/mnt/spaldingdata/nate/communications/coldGrant/BWmog/');
    myClickTracker('/mnt/spaldingdata/nate/communications/coldGrant/BWB73/');
    myClickTracker('/mnt/spaldingdata/nate/communications/coldGrant/BWMo17/');
    myClickTracker('/mnt/spaldingdata/nate/mirror_images/maizeData/sstelpflug/root hair pilot images');
%}