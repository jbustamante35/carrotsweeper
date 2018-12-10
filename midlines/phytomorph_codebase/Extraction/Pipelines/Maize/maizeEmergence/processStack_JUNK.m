function [] = processEmergenceStack(imageStack)

    N = 50;
    rec = getRectification(imageStack{1});
    circleMatrix = getCircles(imageStack,rec,N);
    
    
end

%{

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) scan for images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FilePath = '/home/nate/Downloads/Overhead_Compilation/';
    FilePath = '/home/nate/Downloads/Angle_Compilation/';
    FilePath = '/home/nate/Downloads/emergance/';
    FilePath = '/home/nate/Downloads/20151222_Camera1/';
    FileList = {};
    FileExt = {'tiff','TIF','tif','JPG','jpg'};
    verbose = 1;
    FileList = gdig(FilePath,FileList,FileExt,verbose);
%}