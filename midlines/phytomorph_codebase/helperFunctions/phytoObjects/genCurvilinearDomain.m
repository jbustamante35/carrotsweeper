function [Domain] = genCurvilinearDomain(iCurve,frameRho,WIDTH,nP,I)
    
    %%%%%%%%%%%%%%%%%%%%%%
    % get ref frames along midline    
    BV = igetFrame(iCurve,frameRho);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % create physioframe
    para{1}.type = 'box';
    para{1}.value{1} = [0 0  1];
    para{1}.value{2} = [-WIDTH WIDTH nP];    
    para = genDomains(para);
    
   

    disp = 0;
    Domain = zeros(size(iCurve,1),2,nP);
    
    for pt = 1:size(iCurve,1)
        %%%%%%%%%%%%%%%%%%
        %fprintf(['Start@' num2str(pt) ':' num2str(size(iCurve,1)) '\n']);
        %%%%%%%%%%%%%%%%%%
        % create affine transform
        curPoint = iCurve(pt,:);                        
        %T = [imrotate(squeeze(BV(pt,:,:)),90) fliplr(curPoint)'];
        T = [squeeze(BV(pt,:,:))' fliplr(curPoint)'];
        Domain(pt,:,:) = T*para{1}.d;
        %%%%%%%%%%%%%%%%%%
        % disp
        if disp
            qLEN = 10;
            imshow(I,[]);
            hold on
            plot(iCurve(:,2),iCurve(:,1),'b.')
            plot(iCurve(:,2),iCurve(:,1),'g');
            quiver(iCurve(:,2),iCurve(:,1),qLEN*BV(:,1,2),qLEN*BV(:,2,2),0,'Color','b');
            quiver(iCurve(:,2),iCurve(:,1),qLEN*-BV(:,1,2),qLEN*-BV(:,2,2),0,'Color','r');
            quiver(iCurve(pt,2),iCurve(pt,1),qLEN*-BV(pt,1,2),qLEN*-BV(pt,2,2),0,'Color','g');
            hold off
            drawnow
        end
        %%%%%%%%%%%%%%%%%%
        %fprintf(['End@' num2str(pt) ':' num2str(size(iCurve,1)) '\n']);
    end
    % Domain : [curve x normal] --> R2
    Domain = permute(Domain,[1 3 2]);

end
            