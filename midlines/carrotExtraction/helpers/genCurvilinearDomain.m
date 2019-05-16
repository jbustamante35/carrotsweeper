function domain = genCurvilinearDomain(crv, rho, wid, np, msk, vis)
%% genCurvilinearDomain:
% This function does something and uses igetFrame and genDomains, which uses the
% sub-function createDomains.
%
% Usage:
%   domain = genCurvilinearDomain(crv, RHO, WIDTH, nP, img, vis)
%
% Input:
%   crv: x-/y-coordinates
%   rho: constant representing some value
%   wid: is this a constant too?
%   np: I don't even know
%   img: grayscale image corresponding to curve
%   vis: boolean to visualize output
%
% Output:
%   domain: it's not clear what this is but it's the output
%

%% Get ref frames along midline
BV = igetFrame(crv,rho);

%% Create physioframe
para{1}.type = 'box';
para{1}.value{1} = [0 0  1];
para{1}.value{2} = [-wid wid np];
para = genDomains(para);

domain = zeros(size(crv,1),2,np);

%%
for pt = 1:size(crv,1)
    
    %fprintf(['Start@' num2str(pt) ':' num2str(size(iCurve,1)) '\n']);
    
    %% create affine transform
    curPoint = crv(pt,:);
    %T = [imrotate(squeeze(BV(pt,:,:)),90) fliplr(curPoint)'];
    T              = [squeeze(BV(pt,:,:))' fliplr(curPoint)'];
    domain(pt,:,:) = T * para{1}.d;
    
    %% Visualize output
    if vis && mod(pt/10,1) == 0
        qLEN = 20;
        %         imshow(msk, []);
        %         hold on;
        %         plot(crv(:,2),crv(:,1),'b.');
        %         plot(crv(:,2),crv(:,1),'g');
        %         quiver(crv(pt,2), crv(pt,1), qLEN * -BV(pt,1,2), qLEN * -BV(pt,2,2), ...
        %             0, 'Color', 'g');
        %         waitforbuttonpress;
        
        imagesc(msk);
        axis auto;
        hold on;
        plt(crv, 'm.', 3);
        px = crv(pt, 1);
        py = crv(pt, 2);
        bx = qLEN * -BV(pt, 2, 2);
        by = qLEN * -BV(pt, 1, 2);
        quiver(px, py, bx, by, 0, 'Color', 'r');
        quiver(px, py, by, bx, 0, 'Color', 'b');
        hold off;
        pause(0.001);        
    end
    
    %fprintf(['End@' num2str(pt) ':' num2str(size(iCurve,1)) '\n']);
end

%% Domain : [curve x normal] --> R2
domain = permute(domain, [1 3 2]);

end
