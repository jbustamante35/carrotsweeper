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
    %{
    if (sum(s(:)) == 1) & (s(5) == 0)
        P(e) = 1;
    end
    %}
end


MTI = cat(2,logical(T),~logical(T));
fidx = find(P);
%%
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
tmp = imcrop(tmp,[]);

%{
tmp = [[0 1 0];[1 0 1];[0 1 0]];
tmp = padarray(tmp, [30 20], 0, 'both');
%}

SZ = size(tmp);
tmp = padarray(tmp, [1 1], 'replicate', 'both');


p0 = interp1(xi0,f0,tmp,'linear',eps);
p1 = interp1(xi1,f1,tmp,'linear',eps);



%{
p1 = .8*tmp;
%p0 = double(~logical(p1));
p0 = imcomplement(p1);
%}

T = p0 + p1;
p1 = p1 .* T.^-1;
p0 = p0 .* T.^-1;


p_O_0 = p0;
p_O_1 = p1;


%p_O_1 = imcomplement(p_O_0);


midx0 = P == 0;
midx1 = P == 1;
close all
PS = [];



recent = 0;


for e = 1:2000000
    
    
    [snp0 snp1] = update(p0,p1,MTI,midx0,midx1,SZ,p_O_0,p_O_1);
    
    
    
    imshow(p1,[]);
    drawnow
   
    
    e
    
    
    
    deltaP = abs(p1 - snp1);
    
    pidx(e) = recent(1);
    while sum(recent == pidx(e)) ~= 0
        [V(e) pidx(e)] = max(deltaP(:));
        deltaP(pidx(e)) = 0;
    end
    
    
    p1(pidx(e)) = snp1(pidx(e));
    p0(pidx(e)) = snp0(pidx(e));
    %{
    MD  = (snp1 - p1);
    
    
    
    c1 = p1;
    c0 = p0;
    Z = zeros(size(p1));
    Z(pidx(e)) = 1;
    
    fidx = 1;
    while ~isempty(fidx)
        %z = 1:100
        tZ = imdilate(Z,strel('disk',z,0));
        
        [tsnp0 tsnp1] = update(c0,c1,MTI,midx0,midx1,SZ,p_O_0,p_O_1);
        
        mD = (tsnp1 - c1);
        
        fD = (mD - MD);
        
        
        fidx = find(tZ.*fD ~= 0);
        
        
        %c0(fidx) = tsnp0(fidx);
        %c1(fidx) = tsnp1(fidx);
        
        c0(fidx) = c0(fidx) - fD(fidx);
        c1(fidx) = c1(fidx) + fD(fidx);
        if ~isempty(fidx)
            imshow(c1,[]);
            %waitforbuttonpress
            drawnow;
        end
    end
    
    %}
    
    Z = zeros(size(p1));
    Z(pidx(e)) = 1;
    Z = imdilate(Z,ones(11));
    out = flattenMaskOverlay(p1,logical(Z),1);
    imshow(out,[]);
    drawnow
    %waitforbuttonpress
    
    recent = [recent pidx(e)];
    if numel(recent) > 100
        recent(1) = [];
    end
    pidx(e)
end
%%
close all
for e = 1:size(PS,3)
    imshow(PS(:,:,e),[]);
    drawnow
    waitforbuttonpress
end
%%
close all
for e = 1:size(PS,3)
    imshow(PS(:,:,e),[])
    drawnow
    pause(.5)
end




























