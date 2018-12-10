function [] = specialLoader()
    %%%%%%%%%%%%%%%%%%%%
    % special loader for mat datasets and corresponding thumb    
    %%%%%%%%%%
    % find thumb
    %%%%%%%%%%
    imgPath = '/mnt/scratch5/dev_doane_run2/thumb/';
    imgList = {};
    fileExt = {'jpg','JPG','tif','TIFF','tif'};
    verbose = 1;
    imgList = gdig(imgPath,imgList,fileExt,verbose);
    
    %%%%%%%%%%
    % find patch
    %%%%%%%%%%
    imgPath = '/mnt/scratch5/dev_doane_run2/patch/';
    matList = {};
    fileExt = {'mat'};
    verbose = 1;
    matList = gdig(imgPath,matList,fileExt,verbose);
    
    vwPatch = 0;
    selPts  = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over selection
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    model = 'null';      
    data.domain = [];
    data.codomain = [];
    
    totN = numel(matList);
    totN = 30;
    
    for e = 1:totN
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data
        d = load(matList{e});
        % load thumb
        I = imread(imgList{e});
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % view patch
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if vwPatch
            for p = 1:size(d.O.S,ndims(d.O.S))
                tmp = d.O.S(:,:,:,p);
                imshow(tmp,[]);
                title(num2str(p));
                drawnow
                pause(1);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape data an stack domain and codomain
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % color to gray convert vector
        M = [0.2989 0.5870 0.1140];     
        t.d = double(d.O.S)/255;
        t.s = size(t.d);
        t = nmmult(t,M,3);
        t.d = squeeze(t.d);
        t.s = size(t.d);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % fold time and sample together
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        t.d = reshape(t.d,[t.s(1:2) prod(t.s(3:end))]);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % run patchBalancedPCA over patches (pbPCA)
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate domain
        clear para;
        para{1}.type = 'disk';
        para{1}.value{1} = [0 45 45];
        para{1}.value{2} = [-pi pi 100];
        para = genDomains(para);
        
        
        new_patchSet = pbPCA_loop(t.d,para{1}.d);        
        nP = reshape(new_patchSet,[para{1}.sz t.s(3:4)]);
        
        
        tmp.domain = new_patchSet';
        tmp.domain = hilbertNormalize(tmp.domain);
        tmp.codomain = [];
        
        %idx = data.codomain==1;
        %G = reshape(data.domain(idx,:)',[para{1}.sz sum(idx)]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % gen hypothesis
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        [model grade] = initEvalModel(tmp,model);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % selector
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if selPts
            tP = d.O.para.thumb.percent.value*d.O.P;            
            [pS selV] = mySelector(I,tP,tP(grade,:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % stack data
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        data.domain = [data.domain;tmp.domain];
        data.codomain = [data.codomain;selV];
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % regenerate model
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        model = 'null';        
        [model grade] = initEvalModel(data,model);
    end
end