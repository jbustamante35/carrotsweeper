function [] = tracePath(fileList,PATHS)
    I = imread(fileList{1});
    for sp = 1:numel(PATHS)
        
        
        for pth = 1:numel(PATHS{sp}{1})
            initCurve = PATHS{sp}{1}{pth};
            imshow(I,[]);
            hold on
            for tm = 2:numel(PATHS{sp})
                
                
                d = [];
                for eachP = 1:numel(PATHS{sp}{tm})
                    testPath = PATHS{sp}{tm}{eachP};
                    
                    d(eachP) = pathDistance(initCurve,testPath);
                end
                
                
                
                
                [~,selP(tm)] = min(d);
                
                
                nextCurve = PATHS{sp}{tm}{selP(tm)};
                imshow(I,[])
                hold on
                plot(initCurve(:,2),initCurve(:,1))
                plot(nextCurve(:,2),nextCurve(:,1))
                hold off
                waitforbuttonpress
                
                
                
                initCurve = nextCurve;
            end
            
            
            imshow(I,[]);
            hold on
            plot(initCurve(:,2),initCurve(:,1),'b')
            for tm = 2:numel(PATHS{sp})
                foundP = PATHS{sp}{tm}{selP(tm)};
                plot(foundP(:,2),foundP(:,1),'r')
            end
            waitforbuttonpress
            hold off
            
        end
    end
end