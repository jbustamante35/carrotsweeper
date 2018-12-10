function [bkSIM] = generateBackground(vB1,hB1,vB2,hB2,sz)
    for k = 1:3
        nvB1(1,:,k) = interp1(1:size(vB1,2),squeeze(vB1(:,:,k)),linspace(1,size(vB1,2),sz(2)));
        nvB2(1,:,k) = interp1(1:size(vB2,2),squeeze(vB2(:,:,k)),linspace(1,size(vB2,2),sz(2)));
    end
    
    for k = 1:3
        nhB1(:,1,k) = interp1(1:size(hB1,1),hB1(:,1,k),linspace(1,size(hB1,1),sz(1)));
        nhB2(:,1,k) = interp1(1:size(hB2,1),hB2(:,1,k),linspace(1,size(hB2,1),sz(1)));
    end

    nvB1 = repmat(nvB1,[sz(1) 1 1]);
    nhB1 = repmat(nhB1,[1 sz(2) 1]);
    nvB2 = repmat(nvB2,[sz(1) 1 1]);
    nhB2 = repmat(nhB2,[1 sz(2) 1]);
    bkSIM1 = nvB1 + nhB1;
    bkSIM2 = nvB2 + nhB2;
    bkSIM = .5*(bkSIM1 + bkSIM2);
end