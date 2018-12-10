function [lam] = mulEigs(M)
    F = (M(:,:,1) + M(:,:,4));
    R = ((M(:,:,1) + M(:,:,4)).^2 - 4*(M(:,:,1).*M(:,:,4) - M(:,:,2).^2));
    lam = cat(3,.5*(F+R),.5*(F-R));
    %{
    sz = size(M);
    M = reshape(M,[prod(sz(1:2)) 2 2]);
    for i =1:size(M,1)
        lam(i,:) = eigs(squeeze(M(i,:,:)));
    end
    %}
end