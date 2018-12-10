function [cropBoxes] = applyStream_ver0(FileList,T,num,PARA)
    
    [returnI,cropBoxes] = cropImages_v2(FileList,T{num}{1}.function{1},T{num}{1}.function{2},PARA);
    I = imread(FileList{1});
    imshow(I,[]);
    drawnow
    hold on
    for e = 1:numel(cropBoxes)
        rectangle('Position',cropBoxes{e},'EdgeColor','r');
        drawnow
    end
    
end

%{
    FileList{1} = '/iplant/home/hirsc213/maizeData/seedlingData/04-May-2017/{Plot_883}{Experiment_44}{Planted_4-24-2017}{SeedSource_DI 2113-2}{SeedYear_2016}{Genotype_B73}{Treatment_Control}{PictureDay_10}.nef';
    cropBoxes = applyStream_ver0(FileList,T,1,{950 250});
%}