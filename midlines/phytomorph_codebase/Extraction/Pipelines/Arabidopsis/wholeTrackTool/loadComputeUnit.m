function [UNIT] = loadComputeUnit(toLoad,SZ)
    d = cFlowLoader(toLoad);
    if ~isempty(d)
        vel1 = reshape(squeeze(d(:,1,:)),[SZ size(d,3)]);
        vel2 = reshape(squeeze(d(:,2,:)),[SZ size(d,3)]);
        vel1 = reshape(vel1,[size(vel1,1) size(vel1,2) 1 size(vel1,3)]);
        vel2 = reshape(vel2,[size(vel2,1) size(vel2,2) 1 size(vel2,3)]);
        vel1 = permute(vel1,[2 1 3 4]);
        vel2 = permute(vel2,[2 1 3 4]);
        %{
        K1 = squeeze(vel2(:,:,:,end));
        K2 = squeeze(vel1(:,:,:,end));
        %imshow(I,[]);hold on;drawnow;
        %plot(K1(:),K2(:),'b.');drawnow;
        K1 = squeeze(vel2(:,:,:,1));
        K2 = squeeze(vel1(:,:,:,1));
        %plot(K1(:),K2(:),'r.');hold off;drawnow;
        %}
    else
        vel1 = zeros(SZ(2),SZ(1),1,80);
        vel2 = zeros(SZ(2),SZ(1),1,80);
    end
    UNIT = cat(3,vel2,vel1);
end