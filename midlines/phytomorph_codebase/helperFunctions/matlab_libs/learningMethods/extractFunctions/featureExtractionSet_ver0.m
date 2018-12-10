function [sD] = featureExtractionSet_ver0(sD,fileList,varargin)

    rawWindowSize = varargin{1};
    kurvatureWindowSize = varargin{2};
    kurvatureScaleValue = varargin{3};
    dilateSZ = varargin{4};
    P = varargin{5};
    
    % create import function : imageFile -> rawFeature Object
    func0 = @(data,e)extractFunction_0(data);
    importFunction = fo(sD,func0,'import',{'rawData'});
    image = fmo(sD,'rawImage','type');
  
    % create sliding window function : rawDate -> slidingWindow
    func1 = @(data,e)extractFunction_1(data,rawWindowSize);
    generateSlidingWindows = fo(sD,func1,'generate_windows',{'sliding_window'});
    
    % create kurvature maps : rawData -> kurmap
    func2 = @(data,e)extractFunction_2(data,kurvatureScaleValue,kurvatureWindowSize);
    kurMaps = fo(sD,func2,'generate_K',{'kurvature_1','kurvature_2'});
    
    % create target maps
    func4 = @(data,e)extractFunction_3(data,e,P,dilateSZ);
    tarMaps = fo(sD,func4,'generate_T',{'target_maps'});
       
    % defines a program not done yet : this will define the i:o types 
    rawData = importFunction(image);
    targetData = tarMaps(rawData);
    slideData = generateSlidingWindows(rawData);
    [kur1 kur2] = kurMaps(rawData);
    
    %% can not share sD over the par for loop - therefore need to have empty object returned and then register after in for loop
    for e = 1:numel(fileList)
        out0{e} = importFunction.run(fileList{e});
        %out0{e}.view();
        %drawnow
        out1{e} = tarMaps.run(out0{e},e);
        %out1{e}.view();
        %drawnow
        out2{e} = generateSlidingWindows.run(out0{e});
        %out2{e}.view();
        %drawnow
        [out3{e} out4{e}] = kurMaps.run(out0{e});
        %out3{e}.view();
        %drawnow
    end
    %{
    %% deposit empty objects into feature bank
    for e = 1:numel(fileList)
        sD.depositFeatureMaps(out0{e});
        sD.depositFeatureMaps(out1{e});
        sD.depositFeatureMaps(out2{e});
        sD.depositFeatureMaps(out3{e});
        sD.depositFeatureMaps(out4{e});
    end
    %}
    %check = 1;
    
end