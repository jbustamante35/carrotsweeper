function [data] = myDraw(varargin)
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msgLoop{1} = 'please click';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assign input vars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [filename, extraArgs, msg] = parse_inputs(varargin);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image
    % if I isa cell then only read the first image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isa(filename,'struct')
        I = char(filename.inPort.stack.get(0).getFullFileName()); 
    elseif isa(filename,'cell')
        I = imread(filename{1});
    elseif isa(filename,'char')    
        %{
        info = imfinfo(filename);        
        ROW = [1 samp info.Height];
        COL = [1 samp info.Width];        
        I = imread(filename,'PixelRegion',{ROW,COL});
        %}
        I = myReader(filename);
    elseif isa(filename,'double') & numel(filename) == 1
        I = [];
    else
        I = filename;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get message box
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(extraArgs,'message')
        msgLoop = extraArgs.message;
    end

    for e = 1:numel(msgLoop)
        h1 = msgbox(msgLoop{e},'phytoG','modal');
        uiwait(h1);
        %close h1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample data from user
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if ~isempty(I)
        %    h = figure;
        %    imshow(I);
        %end
        h = figure;
        imshow(I,[]);
        e
        set(0,'CurrentFigure',h);
        drawnow
        pause(.1)
        [x2,x1,V] = impixel();
        data{e} = [x2 x1];
        %close h;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sample data from user
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(extraArgs,'outType')
        convertData = myHS_X('phytoApoint');
        for e = 1:size(data{1},1)
            convertData{e} = extraArgs.outType();
            convertData{e}.setData([data{1}(e,:) 1]);
        end
        data = convertData;
    end
    

   
end




%%%
%%% Function parse_inputs
%%%
function [filename, extraArgs, msg] = parse_inputs(v)
    try
        filename = '';
        extraArgs = {};    
        msg = '';
        % Parse arguments based on their number.
        switch(numel(v))
            case 0
                % Not allowed.
                msg = 'Too few input arguments.';
                return;
            case 1
                % Filename only.
                filename = v{1};
            otherwise
                % Filename and format or other arguments.
                filename = v{1};
                v(1) = [];

                % loop over other arguments
                for e = 1:(numel(v)/2)
                    prop = (e-1)*2 + 1;
                    value = prop + 1;
                    extraArgs.(v{prop}) = v{value};
                end
        end
    catch
        
    end
end



%{
    
    % dig for data for example call
    filePath = '/mnt/spaldingimages/nate/whole_TaN/';
    fileList = {};
    fileExt = {'tiff','TIF'};
    verbose = 1;
    fileList = gdig(filePath,fileList,fileExt,verbose);    



    nI = 3;
    r = round(1 + (numel(fileList)-1).*rand(nI,1));


    bundleRadius = 7;
    radiusSegment = 15;
    % def sample parameters
    sampleParameter{1}.type = 'disk';
    sampleParameter{1}.value{1} = [0 35 35];
    sampleParameter{1}.value{2} = [-pi pi 100];


    data.domain = [];
    data.codomain = [];
    for i = 1:numel(r)
        fn = fileList{r(i)};
        curve = myDraw(fn);
        curve = arcLength(curve,'arcLen',1);
        curveBundle = igetFrame(curve,bundleRadius);
        patchBundle = curveSample(fn,curve,curveBundle,sampleParameter);
        curveSegment = igetCurveSegment(curve,bundleRadius,radiusSegment);
        patchBundle = patchBundle(1:size(curveSegment,1),:,:);
        data.domain = cat(1,data.domain,patchBundle);
        data.codomain = cat(1,data.codomain,curveSegment);
    end


    clear classes;
    load('/home/nate/Desktop/matlab.mat');
    %%%
    d = myT(permute(data.domain,[2 3 1]));
    sD = d.subset(1:10);
    c = myT(permute(data.codomain,[2 3 1]));
    M = myM(d,c);

    %%% create cluster mapping function
    func_cluster{1} = @(x)kmeans(x,2);
    func_cluster{2} = @(x)kmeans(x,3);
    cM = my_cm(func_cluster);
    %%% create least squares fit function
    func_map = my_lsf();
    %%%
    M.addClusterMethod(cM);
    M.addFunctional(func_map);
    M.cluster();
    M.decompose();
    M.i();
    M.e(sD);
    

    t = myT(permute(data.domain,[2 3 1]));
    [model grade] = initEvalModel(data,model)


%}





    
    %{
    %%%%%%%%%%%%%%%
    % arg1 = image
    % arg2 = downSample
    if nargin >= 1;I=varargin{1};end
    if nargin >= 2;samp=varargin{2};else samp=1;end
    
    %%%%%%%%%%%%%%%
    % parse inputs
    [filename, extraArgs, msg] = parse_inputs(varargin);
    
    
    
    if isfield(extraArgs,'message')
        
    end
     
    
    %%%%%%%%%%%%%%%
    % if I isa cell then only read the first image
    if isa(I,'struct')
        I = char(I.inPort.stack.get(0).getFullFileName());
    end
    
    %%%%%%%%%%%%%%%
    % if I isa charater then read
    if isa(I,'char')    
        info = imfinfo(I);        
        ROW = [1 samp info.Height];
        COL = [1 samp info.Width];        
        I = imread(I,'PixelRegion',{ROW,COL});
    end
    
    %%%%%%%%%%%%%%%
    % sample curve from user
    [x2 x1 V] = impixel(I);
    
    curve = [x2 x1];
%}