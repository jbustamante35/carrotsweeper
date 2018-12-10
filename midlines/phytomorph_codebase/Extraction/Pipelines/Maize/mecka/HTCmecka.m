function [] = HTCmecka(user,algorithm,auth)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                HTmecka.m is create jobs for condor to run Maize Ear Cob Kernel Analysis. A user is 
                to choose one of three algorithm to use and provide user for condor.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                mecka.m, ScanAndIssueNewFilesOniRods.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                user:           The user name for condor. (?)
                algorithm:      The argorithm to use. 'c' for singleCobImage.m, 'e' for singleEarImage.m, and 'k' for singleKernelImage.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}


    %{
    %%% start with databse
    import com.franz.agraph.repository.AGServer
    import org.openrdf.model.vocabulary.RDF
    import org.openrdf.model.vocabulary.OWL
    import org.openrdf.model.vocabulary.RDFS
    import org.openrdf.query.QueryLanguage
    %%% connection data
    SERVER_URL = 'http://localhost:10035';
    CATALOG_ID = 'pipelines';
    REPOSITORY_ID = 'functionCalls';
    USERNAME = 'devUser';
    PASSWORD = 'devUser';
    %%% connect
    server = AGServer(SERVER_URL, USERNAME, PASSWORD);
    server.listCatalogs()
    catalog = server.getRootCatalog();
    myRepository = catalog.createRepository(REPOSITORY_ID);
    %%% init
    myRepository.initialize();
    myRepository.isWritable();
    %%% get connection
    conn = myRepository.getConnection();
    %%% rendered needed classes to database
    functionObjects.renderClassDefLevels('farg',conn);
    farg.renderClass(conn);
    fargList.renderClass(conn);
    %}

    % grant access to phytomorphuser
    ga = ['/iplant/home/' user '/#plantType#/'];
    cmd = ['ichmod -r inherit ' strrep(ga,'#plantType#','maizeData')];
    system(cmd);
    cmd = ['ichmod -r read phytomorphuser ' strrep(ga,'#plantType#','maizeData')];
    system(cmd);
    cmd = ['ichmod -r write phytomorphuser ' strrep(ga,'#plantType#','maizeData')];
    system(cmd);

    tmpFileLocation = '/mnt/scratch1/maizePipeline/mecka/tmpSubmitFiles/';
    remoteOutputLocation = ['/iplant/home/' user '/#plantType#/return/#tissueType#/'];
    remoteOutputLocation = strrep(remoteOutputLocation,'#plantType#','maizeData');
    switch algorithm
        case 'c'
            analysisType = 'cobs';
            memREQ = '4000';
            algorithmFlag = 'c';
            numberOfObjects = '3';
            imageRES = '1200';
            localOutputLocation = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/cobData/'];
            remoteOutputLocation = strrep(remoteOutputLocation,'#tissueType#','cobData');
        case 'e'
            analysisType = 'ears';
            memREQ = '10000';
            algorithmFlag = 'e';
            numberOfObjects = '3';
            imageRES = '1200';
            localOutputLocation = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/earData/'];
            remoteOutputLocation = strrep(remoteOutputLocation,'#tissueType#','earData');
        case 'k'
            analysisType = 'kernels';
            memREQ = '4000';
            algorithmFlag = 'k';
            numberOfObjects = [];
            imageRES = '1200';
            localOutputLocation = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/kernelData/'];
            remoteOutputLocation = strrep(remoteOutputLocation,'#tissueType#','kernelData');
    end
    CMD = ['imkdir -p ' remoteOutputLocation];
    system(CMD);
    
    % get file list, transform file list for URL, issue tickets for write
    [FileList rawFileList] = ScanAndIssueNewFilesOniRods(user,analysisType,'maize',{'tif','TIF','tiff','nef'},1);
    %[FileList] = xform2URL(FileList);
    numJobs = numel(FileList);
    [remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numJobs,'write');
    
    %numJobs = 10;
    
    
    func = cFlow('mecka',['maize' '-' analysisType],'1.0');
    %func.setMCRversion('v840');
    func.setMCRversion('v930');
    func.setMemory(memREQ);
    
    %{
        for e = 1:numel(FileList)
            rawFileList{e} = strrep(rawFileList{e},'/iplant/home/nhaase/maizeData/','http://Z.Z.Z.Z/');
            [pth,nm,ext] = fileparts(rawFileList{e});
            FileList{e} = [nm ext];
        end
    %}
    
    
    
    %numJobs = 10;
    for e = 1:numJobs
        fprintf(['start generating job:' num2str(e) ':' num2str(numJobs) '\n']);
        %mecka(algorithmFlag,FileList{e},numberOfObjects,'./output/',remoteOutputLocation,1,1,imageRES,1);
        func(rawFileList{e},algorithmFlag,FileList{e},numberOfObjects,'./output/',remoteOutputLocation,1,1,imageRES,1);
        fprintf(['end generating job:' num2str(e) ':' num2str(numJobs) '\n']);
    end
    %{
    if numJobs ~= 0
        [CMD] = uncLog('ph:l:',rawFileList(1:numJobs),'add',['maize' '-' analysisType],'1.0',{'0'},{'1'},1);
    end
    %}
    
    
    
    
    func.submitDag(auth,150,150);
    
    
    
end

%{

    %auth = csvread('/mnt/spaldingdata/nate/auth.iplant');
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    HTCmecka('garf0012','e');
    HTCmecka('garf0012','c');
    HTCmecka('gxe','c',auth);
    HTCmecka('gxe','e');
    HTCmecka('gxe','k');
    HTCmecka('nhaase','e');
    HTCmecka('nhaase','c');
    HTCmecka('nhaase','k');
    HTCmecka('garf0012','e',auth);
    HTCmecka('garf0012','c',auth);
    HTCmecka('garf0012','k',auth);
    HTCmecka('guoshengwu','e');
    HTCmecka('mocklerlab','k');
    HTCmecka('mrwhite4','e',auth);
    HTCmecka('mrwhite4','c',auth);
    HTCmecka('mrwhite4','k',auth);
    HTCmecka('xinxinding92','k',auth);
    HTCmecka('kmichel','c',auth);
    HTCmecka('kmichel','k',auth);
    HTCmecka('kmichel','e',auth);




    mecka(algorithmFlag,FileList{1},numberOfObjects,'./output/',remoteOutputLocation,1,1,imageRES,1);







    % single runs
    FileList = {'/iplant/home/mocklerlab/maizeData/kernelData/GH3_629_HP301_50.tif'};
    FileList = {'/iplant/home/mrwhite4/maizeData/earData/Scan1-170118-0001.tif'};
   
    mecka('e',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);
    FileList = {'/iplant/home/mocklerlab/maizeData/kernelData/0001582_PI542718.tif'};
    mecka('k',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);

    
    FileList = {'/iplant/home/mocklerlab/maizeData/kernelData/0001582_PI542718.tif'};
    FileList = {'/iplant/home/kmichel/maizeData/cobData/17-06-02/WISN16_110_imageData_cob.tif'};
    mecka('c',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);
   

    FileList = {'/home/nate/Downloads/stevenrbecker/BH_test_cob_01.tiff'};
    mecka('c',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);


    FileList = {'/home/nate/Downloads/husain_agha/TB17_659_b102xb055_cob_test.tiff'};
    mecka('c',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);

    FileList = {'/home/nate/Downloads/stevenrbecker/E_TS_01.tiff'};
    mecka('e',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);

    FileList = {'/iplant/home/mmazaheri/maizeData/earData/17-12-19/BB3-6-1_10Rep1_imageData_ear1.tif'};
    mecka('e',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);

    FileList = {'/iplant/home/mmazaheri/maizeData/kernelData/17-12-19/BB3-6-1-10Rep2_imageData_kernel.tif'};
    mecka('k',FileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);










%}