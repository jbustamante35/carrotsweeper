function domain = genCurvilinearDomain(mline, rho, np, wid, msk, vis)
%% genCurvilinearDomain:
% This function does something and uses igetFrame and genDomains, which uses the
% sub-function createDomains.
%
% Usage:
%   domain = genCurvilinearDomain(mline, rho, np, wid, msk, vis)
%
% Input:
%   mline: x-/y-coordinates
%   rho: constant representing some value
%   np: height of final image [size of image + 1]
%   wid: half the size of np
%   msk: grayscale image corresponding to curve
%   vis: boolean to visualize output
%
% Output:
%   domain: it's not clear what this is but it's the output
%

%% Get ref frames along midline
BV = igetFrame(mline,rho);

%% Create physioframe
para{1}.type     = 'box';
para{1}.value{1} = [0 0 1];
para{1}.value{2} = [-wid wid np];
para             = genDomains(para);
domain           = zeros(size(mline, 1), 2, np);

%%
for pt = 1:size(mline,1)
    %% create affine transform
    curPoint = mline(pt,:);
%     T = [imrotate(squeeze(BV(pt,:,:)),90) fliplr(curPoint)'];
        T              = [squeeze(BV(pt,:,:))' fliplr(curPoint)'];
    domain(pt,:,:) = T * para{1}.d;
    
    %% Visualize output
    if vis && pt == 1
        imagesc(msk);
        colormap gray;
        hold on;
        plt(mline, 'r.', 3);
        ttl = sprintf('Curvilinear Domain');
        title(ttl);
    end
    
    if vis && mod(pt / 10, 1) == 0
        qLEN = 100;
        px   = mline(pt, 1);
        py   = mline(pt, 2);
        bx   = qLEN * -BV(pt, 2, 2);
        by   = qLEN * -BV(pt, 1, 2);
        quiver(px, py, bx, by, 0, 'Color', 'm');
        quiver(px, py, by, bx, 0, 'Color', 'g');
        pause(0.001);
    end
    
end

%% Domain : [curve x normal] --> R2
% domain = permute(domain, [1 3 2]);
domain = permute(domain, [3 1 2]);

end
