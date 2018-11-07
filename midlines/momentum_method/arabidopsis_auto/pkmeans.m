function [U E l CS] = pkmeans(curve,point,n)
    % map each point to its nearest curve point
    for p = 1:size(point,1)
        d = bsxfun(@minus,curve,point(p,:));
        d = sum(d.*d,2);
        [J mp(p)] = min(d);
    end
    
    dC = diff(curve,1,1);
    dL = cumsum([0;sum(dC.*dC,2).^.5]);
    
    
    
    % for each curve point
    for e = 1:size(curve,1)-n-1
        fidx = find((mp >= e) & (mp <= e+n));
        T = mean(dC(e:e+n,:),1);
        T = T/norm(T);
        
        tmpP = bsxfun(@minus,point(fidx,:),curve(e+(n-1)/2,:));
        
        %ontoT = (T'*(point(fidx,:)*T')')';
        ontoT = (T'*(tmpP*T')')';
        subS = point(fidx,:) - ontoT;
        [S C U(e,:) E(:,:,e) L ERR LAM] = PCA_FIT_FULL(subS,3);
        l(:,e) = diag(LAM);
        U(e,:) = curve(e+(n-1)/2,:);
    end
    
    %{
    SQ = zeros(1,size(curve,1));
    SQ(1:n) = 1;
    % for each curve point
    for e = 1:size(curve,1)
        Fidx = [];
        fidx = find(SQ);
        for i = 1:numel(fidx)
            Fidx = [Fidx find(mp == fidx(i))];
        end
        
        T = mean(dC(fidx,:),1);
        T = T/norm(T);
        
        tmpP = bsxfun(@minus,point(fidx,:),curve(round(mean(fidx)),:));
        
        %ontoT = (T'*(point(fidx,:)*T')')';
        ontoT = (T'*(tmpP*T')')';
        subS = point(fidx,:) - ontoT;
        [S C U(e,:) E(:,:,e) L ERR LAM] = PCA_FIT_FULL(subS,3);
        l(:,e) = diag(LAM);
        U(e,:) = curve(e+(n-1)/2,:);
        if e > round(n/2) & e < size(curve,1) - round(n/2)
            SQ = circshift(SQ,[0 1]);
        end
    end
    %}

    
    
    % align vectors
    for e = 2:size(E,3)
        if sign(E(:,1,e)'*E(:,1,e-1)) < 0
            E(:,1,e) = -E(:,1,e);
        end

        if sign(E(:,2,e)'*E(:,2,e-1)) < 0
            E(:,2,e) = -E(:,2,e);
        end
    end
    
    
    
   
   CS = zeros(size(point));
   for e = 1:size(curve,1)-n-1
        ei = e+(n-1)/2;
        fidx = find(mp == ei);
        tmpP = point(fidx,:);
        C = PCA_REPROJ(tmpP,E(:,:,e),U(e,:));
        CS(fidx,2:3) = C(:,1:2);
        CS(fidx,1) = dL(ei);
   end
    
    
end