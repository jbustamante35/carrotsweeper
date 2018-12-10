function [] = collateJSONlocal(sourceLocation,targetLocation)
    massDownload(sourceLocation, '.json',targetLocation);
    FileList= gdig(targetLocation,{},{'json'},1);
    JSONcompile(FileList,targetLocation)
end

%{
    sourceLocation = '/iplant/home/kmichel/maizeData/return/cobData/';
    targetLocation = '/home/nate/Downloads/tempKathrynJSON/cobs/';

    sourceLocation = '/iplant/home/kmichel/maizeData/return/earData/';
    targetLocation = '/home/nate/Downloads/tempKathrynJSON/ears/';


    collateJSONlocal(sourceLocation,targetLocation);
%}