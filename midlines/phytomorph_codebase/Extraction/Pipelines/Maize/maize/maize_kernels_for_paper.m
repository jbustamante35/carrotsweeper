    IDX = [];
    for e = 1:numel(HEADERS)
        if ~isempty(strfind(HEADERS{e},'Kernel Width@ Position'))
        IDX = [IDX e];
        end
    end
    vec = cell2mat(data(2:end,IDX));
    l = cell2mat(data(2:end,9077));
    l(rmidx) = [];
    for e = 1:numel(l)
        L = linspace(-l(e)/2,l(e)/2,size(vec,2));
        CON(:,:,e) = [vec(e,:)' L'];
        e
    end
    plot(CON(:,1,1),CON(:,2,1))
    hold on
    plot(-CON(:,1,1),CON(:,2,1))
    sz = size(CON);
    CONr = reshape(CON,[prod(sz(1:end-1)) sz(end)]);
    rm = any(isnan(CONr),1);
    CONr(:,rm) = [];
    [S C U E L ERR LAM] = PCA_FIT_FULL_T(CONr,size(CONr,1));
    
    CL = {'r' 'g' 'b' 'k' 'y' 'm' 'c' 'r' 'g' 'b'};
    close all
    uC = mean(C,2);
    for e = 1:3
        tmpC = uC;
        tV = linspace(min(C(e,:)),max(C(e,:)),3);
        figure;
        H = [];
        for f = 1:numel(tV)
            tmpC(e) = tV(f);
            M = PCA_BKPROJ_T(tmpC,E,U);
            M = reshape(M,sz(1:end-1));
            plot(M(:,1),M(:,2));
            hold on
            plot(-M(:,1),M(:,2));
            H = [H flipud(M(:,1)) M(:,2) flipud(-M(:,1)) M(:,2)];
        end
        csvwrite(['/mnt/spaldingdata/nate/communications/papers/maizeEarScan/data/kernelPCA_' num2str(e) '.csv'],H)
    end
    %%
    IDX1_1 = [];
    IDX2_1 = [];
    IDX3_1 = [];
    IDX1_2 = [];
    IDX2_2 = [];
    IDX3_2 = [];
    for e = 1:numel(HEADERS)
        if ~isempty(strfind(HEADERS{e},'Ear-1 Width@') & strfind(HEADERS{e},'Scan1'))
            IDX1_1 = [IDX1_1 e];
        end
         if ~isempty(strfind(HEADERS{e},'Ear-2 Width@') & strfind(HEADERS{e},'Scan1'))
            IDX2_1 = [IDX2_1 e];
         end
         if ~isempty(strfind(HEADERS{e},'Ear-3 Width@') & strfind(HEADERS{e},'Scan1'))
            IDX3_1 = [IDX3_1 e];
         end
         if ~isempty(strfind(HEADERS{e},'Ear-1 Width@') & strfind(HEADERS{e},'Scan2'))
            IDX1_2 = [IDX1_2 e];
        end
         if ~isempty(strfind(HEADERS{e},'Ear-2 Width@') & strfind(HEADERS{e},'Scan2'))
            IDX2_2 = [IDX2_2 e];
         end
         if ~isempty(strfind(HEADERS{e},'Ear-3 Width@') & strfind(HEADERS{e},'Scan2'))
            IDX3_2 = [IDX3_2 e];
        end
    end
    
    
    
    vec1_1 = (data(2:end,IDX1_1));
    vec1_2 = (data(2:end,IDX1_2));
    
    
    
    vec2_1 = (data(2:end,IDX2_1));
    vec2_2 = (data(2:end,IDX2_2));
    
    
    vec3_1 = (data(2:end,IDX3_1));
    vec3_2 = (data(2:end,IDX3_2));
    
    l1_1 = find(strcmp(HEADERS,'Max Length Ear1-Scan1'));
    l1_2 = find(strcmp(HEADERS,'Max Length Ear1-Scan2'));
    l2_1 = find(strcmp(HEADERS,'Max Length Ear2-Scan1'));
    l2_2 = find(strcmp(HEADERS,'Max Length Ear2-Scan2'));
    l3_1 = find(strcmp(HEADERS,'Max Length Ear3-Scan1'));
    l3_2 = find(strcmp(HEADERS,'Max Length Ear3-Scan2'));
    
    l1_1 = (data(2:end,l1_1));
    l1_2 = (data(2:end,l1_2));
    
    
    
    l2_1 = (data(2:end,l2_1));
    l2_2 = (data(2:end,l2_2));
    
    
    l3_1 = (data(2:end,l3_1));
    l3_2 = (data(2:end,l3_2));
    
    vecT = [];
    vecL = [];
    M1 = zeros(size(vec1_1,1),1000);
    M2 = zeros(size(vec1_1,1),1000);
    M3 = zeros(size(vec1_1,1),1000);
    for e = 1:size(vec1_1)
        try
            tmp1 = cell2mat(vec1_1(e,:));
            tmp2 = cell2mat(vec1_2(e,:));
            l1tmp = cell2mat(l1_1(e));
            l2tmp = cell2mat(l1_2(e));
            utmp = nanmean([tmp1;tmp2]);
            ul = nanmean([l1tmp;l2tmp]);
            %vecT = [vecT;utmp];
            M1(e,:) = utmp;
            vecL1(e) = [ul];
       
        catch
            
        end

           try
            tmp1 = cell2mat(vec2_1(e,:));
            tmp2 = cell2mat(vec2_2(e,:));
            l1tmp = cell2mat(l2_1(e));
            l2tmp = cell2mat(l2_2(e));
            utmp = nanmean([tmp1;tmp2]);
            ul = nanmean([l1tmp;l2tmp]);
            %vecT = [vecT;utmp];
            M2(e,:) = utmp;
            vecL2(e) = [ul];
       
        catch
            
        end

        
       try
            tmp1 = cell2mat(vec3_1(e,:));
            tmp2 = cell2mat(vec3_2(e,:));
            l1tmp = cell2mat(l3_1(e));
            l2tmp = cell2mat(l3_2(e));
            utmp = nanmean([tmp1;tmp2]);
            ul = nanmean([l1tmp;l2tmp]);
            %vecT = [vecT;utmp];
            M3(e,:) = utmp;
            vecL3(e) = [ul];
       
        catch
            
       end
        
       size(vec1_1)
       e
        
    end
    
    vecT = [M1;M2;M3];
    vecL = [vecL1';vecL2';vecL3'];
    rm = any(vecT==0,2) | any(isnan(vecT),2);
    vecT(rm,:) = [];
    vecL(rm) = [];
    for e = 1:numel(vecL)
        L = linspace(-vecL(e)/2,vecL(e)/2,size(vecT,2));
        CON_E(:,:,e) = [vecT(e,:)' L'];
        e
    end
    
    sz = size(CON_E);
    CON_Er = reshape(CON_E,[prod(sz(1:end-1)) sz(end)]);
    rm = any(isnan(CON_Er),1);
    CON_Er(:,rm) = [];
    [Se Ce Ue Ee Le ERRe LAMe] = PCA_FIT_FULL_T(CON_Er,5);
    %%
    toRM = [1:3];
    toRM_N = 50;
    RM = [];
    for e = 1:numel(toRM)
        [J,sidx] = sort(Ce(toRM(e),:));
        RM = [RM sidx(1:toRM_N)];
        RM = [RM sidx(end-toRM_N+1:end)];
    end
    CON_Er(:,RM) = [];
   
    %%
    
    
    CL = {'r' 'g' 'b' 'k' 'y' 'm' 'c' 'r' 'g' 'b'};
    close all
    uCe = mean(Ce,2);
    for e = 1:size(Ce,1)
        tmpC = uCe;
        tV = linspace(min(Ce(e,:)),max(Ce(e,:)),3);
        figure;
        H = [];
        for f = 1:numel(tV)
            tmpC(e) = tV(f);
            M = PCA_BKPROJ_T(tmpC,Ee,Ue);
            M = reshape(M,sz(1:end-1));
            plot(M(:,1),M(:,2));
            hold on
            plot(-M(:,1),M(:,2));
            axis equal
            H = [H flipud(M(:,1)) M(:,2) flipud(-M(:,1)) M(:,2)];
        end
        csvwrite(['/mnt/spaldingdata/nate/communications/papers/maizeEarScan/data/earPCA_' num2str(e) '.csv'],H)
    end
%%
   
    [S C U E L ERR LAM] = PCA_FIT_FULL(vec,3);
    CL = {'r' 'g' 'b' 'k' 'y' 'm' 'c' 'r' 'g' 'b'};
    close all
    uC = mean(C,1);
    for e = 1:3
        tmpC = uC;
        tV = linspace(min(C(:,e)),max(U(:,e)),10);
        figure;
        for f = 1:numel(tV)
            tmpC(e) = tV(f);
            M = PCA_BKPROJ(tmpC,E,U);
            plot(M,1:numel(M),CL{f})
            hold on
            plot(-M,1:numel(M),CL{f})
            hold all;
        end
    end