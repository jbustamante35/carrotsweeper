function [baseFilter f] = generateFilterSet_v2(T,P)
    
    pidx1 = find(P);
    pidx = 1:numel(P);
    %ptms = permn(pidx,3);
    %bidx = reshape(1:81,[9 9]);
    
    % make unique number map
    uqMap = (2*ones(1,9)).^fliplr((0:8));
    
    % look up table for baseCase -> newCase
    % (newCase,baseCase)
    baseFilter = zeros(2^9,numel(pidx));
    
    % look up table for compatability restrictions
    % (restrictions,baseCase,newCase)
    %f1 = zeros(2^9,numel(pidx),2^9);
    %f2 = f1;
    %N = max(sum(T(find(P),:),2))-1;
    N = 2;
    for n = 1:N
        f{n} = zeros(2^9,numel(pidx),2^9);
    end
    
    
    % make index for pot_mask
    idx_val = reshape(1:81,[9 9]);
    % make potential mask
    pot_mask = padarray(ones(5),[2 2],0,'both');
    % make center square mask
    cen_mask = padarray(ones(3),[3 3],0,'both');
    
    tidx = idx_val(find(pot_mask));
    cidx = idx_val(find(cen_mask));
    % loop over each base case
    for b = 1:numel(pidx1)
        % get base case
        baseCase = reshape(T(pidx1(b),:),[3 3]);
        % get the number for the baseCase
        bcN = uqMap*baseCase(:)+1;
        % make mask
        mag_mask = mag(baseCase,3);
        % make realized mask
        real_mask = pot_mask.*mag_mask;
        ridx = find(real_mask);
        % make potential for 5x5 inner grid
        ptms = permn([0 1],numel(ridx));
        
        
        [q1 q2] = find(baseCase);
        kidx = find(baseCase);        
        q1 = q1 - 2;
        q2 = q2 - 2;
        Q = q1.^2 + q2.^2;
        [Q qidx] = sort(Q);
        kidx = kidx(qidx(1:3));
        filterCase = zeros(3);
        filterCase(kidx) = 1;
        
        
        
        %filterCase = baseCase;
        filterCase(5) = 0;
        
        
        
        fidx = find(filterCase);
        
        
        
        
        
        zc = zeros(3);
        zc(5) = 1;
        zc = mag(zc,3);
        zc = imdilate(zc,ones(3));
        z = {};
        zidx = {};
        zidxS = {};
        for e = 1:numel(fidx)
            z{e} = zeros(3);
            z{e}(fidx(e)) = 1;    
            z{e} = mag(z{e},3);
            suby = z{e};
            suby(find(suby)) = 1:9;
            z{e} = z{e}.*zc;
            zidx{e} = idx_val(find(z{e}));
            zidxS{e} = suby(find(z{e}));
        end
        
        %{
        z1 = zeros(3);
        z1(fidx(1)) = 1;
        z2 = zeros(3);
        z2(fidx(2)) = 1;
        
        z1 = mag(z1,3);
        z2 = mag(z2,3);
        
        zidx1 = idx_val(find(z1));
        zidx2 = idx_val(find(z2));
        %}
        
        
        %tempy = repmat(baseFilter(:,bcN),[1 size(ptms,1)
        for c = 1:size(ptms,1)
            tempyI{c} = [];
            tIDX{c} = {};
        end
        tic
        
        parfor c = 1:size(ptms,1)
            
            
            %tmpBF
            tmp = zeros(9);
            tmp(ridx) = ptms(c,:);
            nidx = uqMap*tmp(cidx)+1;            
            %if nidx == 22
                tmpTest = reshape(tmp(tidx),[5 5]);
                %tmpTest_p = padarray(tmpTest,[1 1],0,'both');
                cval = [];
                for comp = 1:numel(zidx)
                    cval(comp) = sum(tmp(zidx{comp}));
                end
                failFlag = any((cval > 2)) | sum(cval==1) ~= 2;
                flag = vTest(tmpTest,P) & ~failFlag;
                if flag
             %       tmp
             %       c
                    nidx = uqMap*tmp(cidx)+1;
                    %baseFilter(nidx,bcN) = 1;
                    
                    tempyI{c} = [nidx];
                    
                    
                    
                    for e = 1:numel(z)
                        pattern = tmp(zidx{e})';
                        subSet = T(:,zidxS{e});
                        tmpF = all(bsxfun(@eq,subSet,pattern),2);
                        %tmpF = find(any(T(:,logical(tmp(zidx{e}))),2));
                        tIDX{c}{e}.ind = [bcN,nidx]
                        tIDX{c}{e}.d = find(tmpF);
                        
                        %f{e}(tmpF,bcN,nidx) = 1;
                    end
                    %{
                    tmpf1 = find(any(T(:,logical(tmp(zidx1))),2));
                    f1(tmpf1,bcN,nidx) = 1;

                    tmpf2 = find(any(T(:,logical(tmp(zidx2))),2));
                    f2(tmpf2,bcN,nidx) = 1;
                    %}
                end
            %end
        end
        
        for c = 1:size(ptms,1)
            if ~isempty(tempyI{c})
                baseFilter(tempyI{c},bcN) = 1;
                for e = 1:numel(z)
                    if ~isempty(tIDX{c})
                        for k = 1:numel(tIDX{c})
                            if ~isempty(tIDX{c}{e})
                                for g = 1:numel(tIDX{c}{e}.d)
                                    f{e}(tIDX{c}{e}.d(g),tIDX{c}{e}.ind(1),tIDX{c}{e}.ind(2)) = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        toc
        fprintf(['done with base case ' num2str(b) ':' num2str(numel(pidx1)) '\n']);
    end
    done = 1;
end