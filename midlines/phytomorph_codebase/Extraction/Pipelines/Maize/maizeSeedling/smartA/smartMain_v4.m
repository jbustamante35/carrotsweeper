function [] = smartMain_v4(fileName,oPath,rPath,cluter_Level0,cluster_Level1,cluster_Level2,cluster_Level3,kE,kU)
    fprintf(['************************************ STARTING STAGE ONE-TWO ************************************ \n']);
    fprintf(['*****************************    image crop mask and skeleton ********************************** \n']);
    [I,returnI,boundingBoxes,connectedMASK,MASK,SKELETON,LEN] = cropImages_ver3(fileName,950,cluter_Level0,cluster_Level1,cluster_Level2,cluster_Level3,kE,kU);
    fprintf(['************************************  ENDING STAGE ONE-TWO  ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************STARTING STAGE THREE************************************ \n']);
    fprintf(['************************************ measure phenotypes ************************************ \n']);
    [phenoTypes,measureBlocks,JSONdoc] = measurePhenotypes(MASK,SKELETON,returnI,LEN);
    fprintf(['************************************ ENDING STAGE THREE ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************STARTING STAGE FOUR ************************************ \n']);
    fprintf(['************************************  display and save  ************************************ \n']);
    displayResults(fileName,I*255,returnI,boundingBoxes,MASK,SKELETON,phenoTypes,oPath,rPath,JSONdoc);
    fprintf(['************************************  ENDING STAGE FOUR ************************************ \n']);
    close all
    
    
end

%{
    %%%%%%%%%%%
    % to publish smartMain
    % run this with variables created from main
    %func = @(X)smartMain(X,.25,funcObject1,funcObject2,funcObject3,'./output/',[],[]);
    %pF = partialFunction(func,'maizeSeedlings');
    %pF.publish();

    %%%%%%%%%%%
    % to publish smartMain_v4
    % run this with variables created from hotfix
    %func = @(X)smartMain_v4(X,'./output/',[],GMModels{3},level2,levelA,levelD,kE,kU);
    %pF = partialFunction(func,'maizeSeedlings_HOTFIX');
    %pF.publish();

    %%%%%%%%%%%
    % bring in old T
    websave(['./maizeSeedlings.mat'],'https://de.cyverse.org/dl/d/7B2CBBBA-EA62-4DA6-836E-8093366BD8BE/maizeSeedlings.mat');
    load('./maizeSeedlings.mat');
    w = functions(obj.func);
    T = w.workspace{1}.T;


    %%%%%%%%%%%
    % to publish smartMain_v2
    % run this with variables created from main
    func = @(X)smartMain_v2(X,.25,T,{1000 250},'./output/',[],[]);
    pF = partialFunction(func,'maizeSeedlings');
    pF.publish();


    % func(FileList{1})
    func('/iplant/home/hirsc213/maizeData/seedlingData/01-May-2017/{Plot_883}{Experiment_44}{Planted_4-24-2017}{SeedSource_DI 2113-2}{SeedYear_2016}{Genotype_B73}{Treatment_Control}{PictureDay_7}.nef')
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/27-Jun-2017/{Plot_1089}{Experiment_52}{Planted_6-14-2017}{SeedSource_DI1910-7}{SeedYear_2015}{Genotype_Mo17}{Treatment_Control}{PictureDay_13}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/23-Jun-2017/{Plot_1095}{Experiment_52}{Planted_6-14-2017}{SeedSource_DI1910-7}{SeedYear_2015}{Genotype_Mo17}{Treatment_1 d 2C8C stress}{PictureDay_9}.nef'
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/22-Jul-2017/{Plot_1277}{Experiment_56}{Planted_7-13-2017}{SeedSource_IBM1697-1}{SeedYear_2012}{Genotype_M345}{Treatment_Control}{PictureDay_9}.nef';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/03-May-2017/{Plot_Hirsch_883}{Experiment_44}{Planted_4-24-17}{SeedSource_DI 2113-2}{SeedYear_2016}{Genotype_B73}{Treatment_Control}{PictureDay_9}.nef';
    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/10-May-2017/{Plot_894}{Experiment_44}{Planted_4-24-2017}{SeedSource_DI 2114-3}{SeedYear_2016}{Genotype_Mo17}{Treatment_Control}{PictureDay_16}.nef';
    %fileName = '/iplant/home/hirsc213/maizeData/seedlingData/05-May-2017/{Plot_Hirsch_883}{Experiment_44}{Planted_4-24-17}{SeedSource_DI 2113-2}{SeedYear_2016}{Genotype_B73}{Treatment_Control}{PictureDay_11}.nef';
    func(fileName); 

    fileName = '/iplant/home/hirsc213/maizeData/seedlingData/20-Dec-2017/{Plot_2040}{Experiment_73}{Planted_12-8-2017}{SeedSource_DI2128-3}{SeedYear_2016}{Genotype_LH185}{Treatment_Control}{PictureDay_12}.nef';
    obj.func(fileName);

%}