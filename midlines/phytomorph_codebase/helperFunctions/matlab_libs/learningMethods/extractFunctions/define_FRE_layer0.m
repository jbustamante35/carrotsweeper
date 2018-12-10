function [fl] = define_FRE_layer0(sD,in,varargin)

    % set the window size
    rawWindowSize = varargin{1};
    
    
    % init FREE layer
    fl = FRE();
    %%%%%%%%%%%%%%%%%%%%%
    % set the feature extraction layer
    %%%%%%%%%%%%%%%%%%%%%
    % create sliding window function : rawDate -> slidingWindow
    func1 = @(data,e)extractFunction_1(data,rawWindowSize);
    generateSlidingWindows = fo(sD,func1,'generate_sliding_window',{'sliding_window'});
    patchData = generateSlidingWindows(in);
    % set the feature extraction layer function
    fl.setF(generateSlidingWindows);
    
    %%%%%%%%%%%%%%%%%%%%%
    % set the reprentation layer 
    fl.setR(@(x)obtainBasisInformation_ver0(x));
    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%
    % set the error layer
    %%%%%%%%%%%%%%%%%%%%%
    func2 = @(data,bi)extractFunction_4(data,bi,rawWindowSize);
    generateErrorMap = fo(sD,func2,'generate_sliding_window_error',{'sliding_window_error'});
    errorData = generateErrorMap(patchData);
    % se the error layer function
    fl.setE(generateErrorMap);
    
    
end