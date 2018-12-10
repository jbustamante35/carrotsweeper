% number of points to interpolate the midline with
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % file list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fileList = sdig([main_path{e} filesep],{},{'TIF','tif','png','PNG'},1);
            fileList = fileList{1};
            %{
            for s = 1:numel(fileList)                
                [pth nm ext] = fileparts(fileList{s});
                N(s) = str2num(nm);                              
            end
            [N sidx] = sort(N);
            fileList = fileList(sidx);
            %}
            %{
            for s = 1:numel(fileList)
                for i = 1:numel(fileList{s})
                    [pth nm ext] = fileparts(fileList{s}{e});
                    N(i) = str2num(nm);
                end
                [N sidx] = sort(N);
                fileList{s} = fileList{s}(sidx);
            end
            %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % threshold for logical statement regarding when to "stop" tracking
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %para.THRESH = 10;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %para.disp = 1;
            prompt = {'Time downsample values (fr):',...
                      'Time window values (fr):',...
                      'Space window values (px):',...
                      'Space downsample values (px):'...
                      'Derivative values (px):'....
                      'QC window size (px):'...
                      'Disk tracking radius (px):'};
            dlg_title = 'Please enter kinematic values:';
            default = {'1','1,3,5,7','30','50,100','15,25,30,40,50','30','70'};             % default down sample for frame
            res = inputdlg(prompt,dlg_title,1,default);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % time downsample value
        % number of frames to skip in the analysis
        % the "helps" or is "needed" becuase, 
        % while the tracking "needs" a high frame rate to 
        % find matching patches, BUT over "short" time intervals,
        % the tracking is noisy.
        % note: a paper NEEDS to be written, OR a paper found
        % which address biological rates, how they are intertwined
        % with the spatio-temporal resolution of the "imaging" device.
        % as well, there is a paper (maybe the same one) which needs to 
        % address a heterogenous analysis methods, for example, the method
        % should "differ" in the expansion zone VS the meristem.
        % AND should differ from trial to trial.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            para(e).frameDownsample = str2num(res{1});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of frames - after downsampling
        % to use in the analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            para(e).frameWindow = str2num(res{2});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % smoothing parameter for derivative analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            para(e).spaceWindow = str2num(res{3});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % smoothing parameter for derivative analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            para(e).space_downsample = str2num(res{4});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % smoothing parameter for derivative analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            para(e).derSMOOTH = str2num(res{5});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % smoothing parameter for derivative analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            para(e).anchor_point_window = str2num(res{6});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % radius of disk,thresh,downsample frames,window
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            para(e).diskRadius = str2num(res{7});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % radius of disk,thresh,downsample frames,window
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            para(e).THRESH = 100;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % outpath
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % main outpath
            % outPath{e} = ['/mnt/scratch3/users/monshausenlab/testRun/ver0/'];
            %outPath{e} = ['/mnt/scratch3/users/trieupham10/results/ver1/'];
            outPath{e} = ['/mnt/scratch3/users/monshausenlab/testRun/ver_2013.05.10/'];
            %outPath{e} = ['/mnt/scratch3/users/nmiller/vfor_E/'];
            % special outpath
            firstFile = fileList{1};
            [pth,nm,ext] = fileparts(firstFile);
            fidx = strfind(pth,filesep);
            nm = pth(fidx(end)+1:end);
            outPath{e} = [outPath{e} nm filesep];  
        
        %}    
            