function [] = specialLoader_bulk()
    %%%%%%%%%%%%%%%%%%%%
    % special loader for mat datasets and corresponding thumb    
    %%%%%%%%%%
    % find thumb
    %%%%%%%%%%
    imgPath = '/mnt/scratch5/dev_doane_run1/thumb/';
    imgList = {};
    fileExt = {'jpg','JPG','tif','TIFF','tif'};
    verbose = 1;
    imgList = gdig(imgPath,imgList,fileExt,verbose);
    
    %%%%%%%%%%
    % find patch
    %%%%%%%%%%
    imgPath = '/mnt/scratch5/dev_doane_run1/patch/';
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
    totN = 10;
    data.domain = [];
    data.location = [];
    data.number = [];
        
    
    
    
    
    for e = 1:totN
        fprintf(['Starting@' num2str(e) ':' num2str(totN) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data
        d = load(matList{e});
        % load thumb
        I = imread(imgList{e});
        
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
        % execute pdPCA over patches
        new_patchSet = pbPCA_loop(t.d,para{1}.d);        
        nP = reshape(new_patchSet,[para{1}.sz t.s(3:4)]);
        
        
       
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % store data
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        data.domain = cat(4,data.domain,permute(nP,[1 2 4 3]));
        data.location = cat(3,data.location,permute(d.O.P,[2 3 1]));
        data.number = cat(1,data.number,e*ones(size(d.O.P,1),1));
        
        
        fprintf(['Ending@' num2str(e) ':' num2str(totN) '\n']);
    end
    
    
    
    
    
    data.selection = [];
    for e = 1:totN
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % load image
        %%%%%%%%%%%%%%%%%%%%%%%%%%        
        % load thumb
        I = imread(imgList{e});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % selector
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        fidx = find(data.number==e);
        if selPts
            tP = d.O.op.para.thumb.percent.value*squeeze(data.location(5:6,1,fidx));
            tP = fliplr(tP');
            [pS selV] = mySelector(I,tP,[]);
        end
        
        data.selection = [data.selection;selV];
    end
    
    
    pck = ipca(data.domain(:,:,:,find(data.selection)));
    t = rotT(data.domain);
    cf = pck.bf'*t;
    sim = pck.bf*cf;
    error = sim-t;
    error = sum(error.*error,1).^.5;
    
    fidx = find(data.selection);
    plot3(cf(1,:),cf(2,:),cf(3,:),'k.')
    hold on
    plot3(cf(1,fidx),cf(2,fidx),cf(3,fidx),'go');
    
    lambda = myLDA(cf',data.selection);
    value = cf'*lambda;
    plot3(value,value,value,'k.');
    hold on
    plot3(value(fidx),value(fidx),value(fidx),'go');
    plot3(value(fidx),value(fidx),value(fidx),'g.');
    
   
    tsz = size(data.domain);
    tsz = tsz(1:end-1);
    g = pck.bf*lambda;
    g = reshape(g,tsz);
    
    data.domain
    hello = 5;
end