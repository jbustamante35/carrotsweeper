function [measureBlocks] = smartMain(fileName,reSize,networkObject1,networkObject2,networkObject3,oPath,rPath,toM)
    fprintf(['************************************ STARTING STAGE ONE ************************************ \n']);
    fprintf(['************************************    image crop      ************************************ \n']);
    % crop the plants with the NN
    I = double(imread(fileName));
    [I,angle] = rectifyImage(I/255);
    I = I*255;
    [returnI,boundingBoxes] = cropImages(I,reSize,networkObject1,networkObject2,networkObject3);
    fprintf(['************************************  ENDING STAGE ONE  ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************ STARTING STAGE TWO ************************************ \n']);
    fprintf(['************************************  mask and skeleton ************************************ \n']);
    % get the mask(s) and skeleton(s)
    [connectedMASK,MASK,SKELETON] = getMASKandSKELETON(returnI);
    fprintf(['************************************  ENDING STAGE TWO  ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************STARTING STAGE THREE************************************ \n']);
    fprintf(['************************************ measure phenotypes ************************************ \n']);
    [phenoTypes measureBlocks] = measurePhenotypes(MASK,SKELETON,toM);
    fprintf(['************************************ ENDING STAGE THREE ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************STARTING STAGE FOUR ************************************ \n']);
    fprintf(['************************************  display and save  ************************************ \n']);
    displayResults(fileName,I,returnI,boundingBoxes,MASK,SKELETON,phenoTypes,oPath,rPath);
    fprintf(['************************************  ENDING STAGE FOUR ************************************ \n']);
    close all
end

%{
    % to publish smartMain
    % run this with variables created from main
    func = @(X)smartMain(X,.25,funcObject1,funcObject2,funcObject3,'./output/',[],[]);
    pF = partialFunction(func,'maizeSeedlings');
    pF.publish();
    % func(FileList{1})
    func('/iplant/home/hirsc213/maizeData/seedlingData/01-May-2017/{Plot_883}{Experiment_44}{Planted_4-24-2017}{SeedSource_DI 2113-2}{SeedYear_2016}{Genotype_B73}{Treatment_Control}{PictureDay_7}.nef')
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/27-Jun-2017/{Plot_1089}{Experiment_52}{Planted_6-14-2017}{SeedSource_DI1910-7}{SeedYear_2015}{Genotype_Mo17}{Treatment_Control}{PictureDay_13}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/23-Jun-2017/{Plot_1095}{Experiment_52}{Planted_6-14-2017}{SeedSource_DI1910-7}{SeedYear_2015}{Genotype_Mo17}{Treatment_1 d 2C8C stress}{PictureDay_9}.nef'
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/22-Jul-2017/{Plot_1277}{Experiment_56}{Planted_7-13-2017}{SeedSource_IBM1697-1}{SeedYear_2012}{Genotype_M345}{Treatment_Control}{PictureDay_9}.nef';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/03-May-2017/{Plot_Hirsch_883}{Experiment_44}{Planted_4-24-17}{SeedSource_DI 2113-2}{SeedYear_2016}{Genotype_B73}{Treatment_Control}{PictureDay_9}.nef';
    func(fileName); 
%}