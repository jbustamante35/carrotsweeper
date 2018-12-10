T = permn([0 1],9);
P = zeros(size(T,1),1);
for e = 1:size(T,1)
    s = reshape(T(e,:),[3 3]);
    m = imfilter(s,ones(2));
    m(:,end) = [];
    m(end,:) = [];
    
    if any(m(:) >= 3)
        P(e) = 0;
    else
        if T(e,5) == 0
            P(e) = 0;
        else
            if sum(T(e,:)) ~= 3
                P(e) = 0;
            else
                P(e) = 1;    
            end
            
        end
    end
end


MTI = cat(2,logical(T),~logical(T));
fidx = find(P);
for e = 1:numel(fidx)
    t = T(fidx(e),:);
    t = reshape(t,[3 3]);
    imshow(t,[]);
    t
    waitforbuttonpress
end

%%
pidx1 = find(sE(:)==1);
pidx0 = find(sE(:)==0);
v0 = E(pidx0);
v1 = E(pidx1);
[f0,xi0] = ksdensity(v0(1:10000),linspace(0,50,500));
[f1,xi1] = ksdensity(v1(1:10000),linspace(0,50,500));
%%
midx0 = T(:,5) == 0;
midx1 = T(:,5) == 1;
%%
tmp = E(:,:,1);
tmp = imresize(tmp,.15);
SZ = size(tmp);
tmp = padarray(tmp, [1 1], 'replicate', 'both');
p0 = interp1(xi0,f0,tmp,'linear',eps);
p1 = interp1(xi1,f1,tmp,'linear',eps);
p_O_0 = p0;
p_O_1 = p1;
midx0 = P == 0;
midx1 = P == 1;
close all
PS = [];

for e = 1:200
    p0(p0==0) = eps;
    p1(p1==0) = eps;
    T = p0 + p1;
    p1 = p1 .* T .^-1;
    p0 = p0 .* T .^-1;

    lp0 = log(p0);
    lp1 = log(p1);
    %lp0(isinf(lp0)) = 0;
    %lp1(isinf(lp1)) = 0;

    %lp0 = padarray(lp0, [1 1], 'replicate', 'both');
    %lp1 = padarray(lp1, [1 1], 'replicate', 'both');
    
    rp0 = im2colF(lp0,[3 3],[1 1]);
    rp1 = im2colF(lp1,[3 3],[1 1]);
    newI = mtimesx(double(MTI),[rp1;rp0]);
    newI0 = newI(midx0,:);
    newI1 = newI(midx1,:);
    
    
    %newI0 = mtimesx((P(midx0)+1),'T',newI0);
    %newI1 = mtimesx(P(midx1),'T',newI1);
    %newI0 = mean(newI0,1)*sum(midx0);
    %newI1 = mean(newI1,1)*sum(midx1);
    
    
    newI1 = max(newI1);
    newI0 = max(newI0);
    
    
    R = newI1 - newI0;
    A = 1;
    %R = bindVec(R);
    R = (exp(R)+1).^-1;
    %{
    MN = max(max(newI0),max(newI1));
    newI0 = (newI0 - MN);
    newI1 = (newI1 - MN);
    
    
    MX = min(min(newI0),min(newI1));
    newI0 = -newI0 / MX;
    newI1 = -newI1 / MX;
    %}
    
    v0 = reshape(R,SZ);
    v1 = 1-v0;
    %{
    v1 = reshape(exp(A*newI1),SZ);
    v0 = reshape(exp(A*newI0),SZ);
    %}
    T = v1 + v0;
    v1 = v1 .* T.^-1;
    v0 = v0 .* T.^-1;

    np1 = padarray(v1, [1 1], 'replicate', 'both');
    np0 = padarray(v0, [1 1], 'replicate', 'both');
    
    %p0 = np0.*p_O_0;
    %p1 = np1.*p_O_1;
    
    p0 = np0.^1.*p_O_0;
    p1 = np1.^1.*p_O_1;
    
    %p0 = .2*np0 + .8*p_O_0;
    %p1 = .2*np1 + .8*p_O_1;
    
    %p0 = np0.*p0;
    %p1 = np1.*p1;
    
    
    %p0 = np0;
    %p1 = np1;
    
    
    T = p1 + p0;
    p1 = p1 .* T.^-1;
    p0 = p0 .* T.^-1;
    
    %{
    T = cat(3,p0,p1);
    [J1 p1] = max(T,[],3);
    p1 = p1 - 1;
    [J0 p0] = min(T,[],3);
    p0 = p0 - 1;
    %}
   
    
    imshow(p1,[]);
    drawnow
    PS(:,:,e) = p1;
end
%%
close all
for e = 1:size(PS,3)
    imshow(PS(:,:,e),[]);
    drawnow
end
%%
close all
for e = 1:size(PS,3)
    imshow(PS(:,:,e),[])
    drawnow
    pause(.5)
end




























