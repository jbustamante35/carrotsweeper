function domain = genCurvilinearDomain(crv, rho, np, wid, msk, vis)
%% genCurvilinearDomain:
% This function does something and uses igetFrame and genDomains, which uses the
% sub-function createDomains.
%
% Usage:
%   domain = genCurvilinearDomain(crv, rho, np, wid, msk, vis)
%
% Input:
%   crv: x-/y-coordinates
%   rho: constant representing some value
%   np: height of final image [size of image + 1]
%   wid: half the size of np
%   img: grayscale image corresponding to curve
%   vis: boolean to visualize output
%
% Output:
%   domain: it's not clear what this is but it's the output
%

%% Get ref frames along midline
BV = igetFrame(crv,rho);

%% Create physioframe
para{1}.type     = 'box';
para{1}.value{1} = [0 0 1];
para{1}.value{2} = [-wid wid np];
para             = genDomains(para);
domain           = zeros(size(crv, 1), 2, np);

%%
for pt = 1:size(crv,1)
    %% create affine transform
    curPoint = crv(pt,:);
%     T = [imrotate(squeeze(BV(pt,:,:)),90) fliplr(curPoint)'];
        T              = [squeeze(BV(pt,:,:))' fliplr(curPoint)'];
    domain(pt,:,:) = T * para{1}.d;
    
    %% Visualize output
    if vis && pt == 1
        imagesc(msk);
        colormap gray;
        hold on;
        plt(crv, 'r.', 3);
        ttl = sprintf('Curvilinear Domain');
        title(ttl);
    end
    
    if vis && mod(pt / 10, 1) == 0
        qLEN = 100;
        px   = crv(pt, 1);
        py   = crv(pt, 2);
        bx   = qLEN * -BV(pt, 2, 2);
        by   = qLEN * -BV(pt, 1, 2);
        quiver(px, py, bx, by, 0, 'Color', 'm');
        quiver(px, py, by, bx, 0, 'Color', 'g');
        pause(0.001);
    end
    
end

%% Domain : [curve x normal] --> R2
domain = permute(domain, [1 3 2]);

end
