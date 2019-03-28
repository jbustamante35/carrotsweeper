function domain = genCurvilinearDomain(crv, RHO, WIDTH, nP, img, vis)
%% genCurvilinearDomain:
% This function does something and uses igetFrame and genDomains, which uses the
% sub-function createDomains.
%
% Usage:
%   domain = genCurvilinearDomain(crv, RHO, WIDTH, nP, img, vis)
%
% Input:
%   crv: x-/y-coordinates
%   RHO: constant representing some value
%   WIDTH: is this a constant too?
%   nP: I don't even know
%   img: grayscale image corresponding to curve
%   vis: boolean to visualize output
%
% Output:
%   domain: it's not clear what this is but it's the output
%

%% Get ref frames along midline
BV = igetFrame(crv,RHO);

%% Create physioframe
para{1}.type = 'box';
para{1}.value{1} = [0 0  1];
para{1}.value{2} = [-WIDTH WIDTH nP];
para = genDomains(para);

domain = zeros(size(crv,1),2,nP);

%%
for pt = 1:size(crv,1)
    
    %fprintf(['Start@' num2str(pt) ':' num2str(size(iCurve,1)) '\n']);
    
    %% create affine transform
    curPoint = crv(pt,:);
    %T = [imrotate(squeeze(BV(pt,:,:)),90) fliplr(curPoint)'];
    T = [squeeze(BV(pt,:,:))' fliplr(curPoint)'];
    domain(pt,:,:) = T*para{1}.d;
    
    %% Visualize output
    if vis
        qLEN = 10;
        imshow(img, []);
        hold on;
        plot(crv(:,2),crv(:,1),'b.');
        plot(crv(:,2),crv(:,1),'g');
        %quiver(iCurve(:,2),iCurve(:,1),qLEN*BV(:,1,2),qLEN*BV(:,2,2),0,'Color','b');
        %quiver(iCurve(:,2),iCurve(:,1),qLEN*-BV(:,1,2),qLEN*-BV(:,2,2),0,'Color','r');
        quiver(crv(pt,2), crv(pt,1), qLEN * -BV(pt,1,2), qLEN * -BV(pt,2,2), ...
            0, 'Color', 'g');
        hold off;
        drawnow;
        waitforbuttonpress;
    end
    
    %fprintf(['End@' num2str(pt) ':' num2str(size(iCurve,1)) '\n']);
end

%% Domain : [curve x normal] --> R2
domain = permute(domain, [1 3 2]);

end
